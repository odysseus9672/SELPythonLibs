# Class to construct binary classification trees and facilitate translation
# between the balanced log-odds of taking a branch on the tree's inner nodes
# and the (log of) the probabilities of being on each of the leaf nodes.
# The log odds are "balanced" because when they are all zero all of the leaf
# node probabilites are the same, regardless of tree shape.
#
# Sean E. Lake @ NAOC
# Aug 18, 2021
# 

import numpy as np
import torch as tch
import torch.nn.functional as tchF

class BinClassTree():
    """
    Translate between a multi-class log-odds list and its corresponding 
    probabilities.

    Instantiantes a binary classification tree. Useful for translating 
    between the balanced log-odds of taking any branch of the tree to the log 
    of the probabilities of each class (and back). The prototypical use-case
    for such a classification tree is to sit between a neural net, or other 
    maachine learning classifier, and translate between ``Nclass - 1`` floats 
    that can take any value and ``Nclass`` non-positive floats that represent 
    the natural log of the probability that the example the classifier is 
    evaluating belongs to that class. The log-odds are "balanced" because when 
    they are all zero then all of the probabilities equal (and the converse).

    This addresses two shortcomings of a soft-max based approach. One, if the 
    full soft-max is taken of the outputs when calculating the loss function 
    then loss of significance and vanishing gradients can become a problem. 
    Two, the soft-max of ``Nclass`` numbers produces ``Nclass`` numbers that 
    have ``Nclass-1`` degrees of freedom, meaning some of the network's 
    resources were wasted computing an irrelevant degree of freedom.

    The structure of the tree is specified by a list of integers that give the 
    split points of the class list. The number of leaf nodes gives the number 
    of classes, and that's one more than the number of split-points. The first
    split point in the list gives the number of leaf nodes on the root node's 
    left branch. The remaining split points are divided between child binary 
    classification trees in such a way that the left child gets the first
    set of points sufficient to specify its structure and the right node gets
    the rest.

    Example trees and their splitpoints:
    /\  = [1]

    /\
     /\ = [1, 1]

     /\
    /\  = [2, 1]

      /\
     /  \
    /\  /\  = [2, 1, 1]

    /\
     /\
      /\ = [1, 1, 1]

      /\
     /\
    /\   = [3, 2, 1]


    Leaves on the left tree: ``[1] * (Nclasses - 1)``
    Leaves on the right tree: ``range(1, Nclasses)[::-1]``
    Further balanced tree examples: ``[4, 2, 1, 1, 2, 1, 1]``
                                    ``[8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]``

    Parameters
    ----------
    splitpoints : sequence, list, tuple, etc of integers
        List of integers specifying the number of leaf nodes below the left
        branch of each inner node. Inner nodes are traversed starting at the 
        root node from left to right.

    """
    
    def __init__(self, splitpoints):
        splitpoints = np.asarray(splitpoints)
        Nclass = len(splitpoints) + 1
        self.Nclass = Nclass
        
        if Nclass == 2:
            if splitpoints[0] != 1:
                raise ValueError("Malformed splitpoint at leaf node.")
            
            self.LOoffset = 0.0
            self.children = (None, None)
            self.Nleft = 1
            
        elif Nclass > 2:
            Nleft = splitpoints[0]
            if Nleft >= len(splitpoints) or Nleft <= 0:
                errtxt = "Splitpoint at {:d} does not split a list length {:d}"
                raise ValueError(errtxt.format(Nleft, Nclass))

            self.Nleft = Nleft
            Nright = Nclass - Nleft
            dN = Nleft - Nright
            if 2*abs(dN) <= Nleft:
                self.LOoffset = np.log1p(-float(dN)/float(Nleft))
            else:
                self.LOoffset = np.log(float(Nright) / float(Nleft))

            childsplits = splitpoints[1:]
            if Nleft == 1:
                self.children = (None, ClassTree(childsplits))
                
            elif Nright == 1:
                self.children = (ClassTree(childsplits), None)
                
            else:
                csp = sp - 1

                self.children = (ClassTree(childsplits[:csp]),
                                 ClassTree(childsplits[csp:]))
        else:
            raise ValueError("length fo splitpoints insufficient")

        return None
    

    def ProbsTologOdds(self, Probs):
        """
        Calculate the balanced log-odds from the classification probabilities.

        Given the probabilities of that a set of examples belongs to each 
        classs, returns the balanced log-odds of teaching each branch of the 
        inner nodes. ``Probs`` must be 2-d, and the second dimension must match 
        the number of classes (leaf nodes) in the tree. The log-odds are 
        balanced in that when all of the probabilities for a given example are 
        equal all of its inner node log-odds are zero.

        Parameters
        ----------
        Probs : `torch.tensor`
            A 2-dimensional `torch.tensor` of classification probabilities. 
            The shape of ``Probs`` must be ``(Nexamples, Nclass)``, where 
            ``Nclass`` is the number of leaf nodes in the binary 
            classification tree.

        Returns
        -------
        logodds : `torch.tensor`
            A 2-dimensional array of balanced inner node log-odds. Will
            have shape ``(Nexamples, Nclass - 1)``.
        """
        if Probs.shape[1] != self.Nclass:
            errtxt = ("Number of probabilities {:d} mismatch number of" +
                      " classes {:d}")
            raise ValueError(errtxt.format(Probs.shape[1], self.Nclass))
        
        Nprobs = Probs.shape[0]
        Nleft = self.Nleft
        results = tch.empty((Nprobs, self.Nclass - 1))
        lo = tch.log(tch.sum(Probs[:,:Nleft], 1) /
                     tch.sum(Probs[:,Nleft:], 1)) + self.LOoffset
        
        results[:,0] = lo
        Rightstart = 1
        if self.children[0] is not None:
            leftprobs = Probs[:,:Nleft]
            results[:,1:self.Nleft] = self.children[0].ProbsTologOdds(leftprobs)
            Rightstart += self.Nleft

        if self.children[1] is not None:
            rightprobs = Probs[:,Nleft:]
            results[:,Rightstart:] = self.children[1].ProbsTologOdds(rightprobs)
        
        return results
    

    def logOddsTologProbs(self, logodds):
        """
        Calculate the log of the classification probabilities from the balanced
        log-odds of the tree's inner nodes.

        Given the balanced log-odds of that a set of examples belongs to each 
        branch of the binary classification tree's inner nodes, returns the 
        log of the probabilities that each example belongs to each class. 
        ``logodds`` must be 2-d, and the second dimension must be one less than
        the number of classes (leaf nodes) in the tree. The log-odds are 
        balanced in that when all of the probabilities for a given example are 
        equal all of its inner node log-odds are zero.

        Parameters
        ----------
        logodds : `torch.tensor`
            A 2-dimensional `torch.tensor` of balanced inner node log-odds. 
            The shape of ``logodds`` must be ``(Nexamples, Nclass - 1)``, 
            where ``Nclass`` is the number of leaf nodes in the binary 
            classification tree.

        Returns
        -------
        probs : `torch.tensor`
            A 2-dimensional `torch.tensor` of the logarithm of classification 
            probabilities. Will have shape ``(Nexamples, Nclass)``.
        """
        if logodds.shape[1] != self.Nclass - 1:
            errtxt = ("Number of logodds {:d} mismatch number of " +
                      "classes {:d} - 1")
            raise ValueError(errtxt.format(logodds.shape[1], self.Nclass))

        Nprobs = logodds.shape[0]
        Nleft = self.Nleft
        results = tch.zeros((Nprobs, self.Nclass))

        netlo = tch.reshape(logodds[:,0] - self.LOoffset, (Nprobs, 1))
        results[:,:Nleft] += tchF.logsigmoid(netlo)
        results[:,Nleft:] += tchF.logsigmoid(-netlo)

        rightstart = 1
        if self.children[0] is not None:
            rightstart = self.children[0].Nleft + 1
            leftlogodds = logodds[:,1:rightstart]
            lchild = self.children[0]
            results[:,:Nleft] += lchild.logOddsTologProbs(leftlogodds)

        if self.children[1] is not None:
            rightlogodds = logodds[:,rightstart:]
            rchild = self.children[1]
            results[:,Nleft:] += rchild.logOddsTologProbs(rightlogodds)

        return results


# Example code:
# SplitHierarchy = (1, 1, 1, 1)
# translator = ClassTree(SplitHierarchy)
