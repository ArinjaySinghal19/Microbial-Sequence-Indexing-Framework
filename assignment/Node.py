#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 23:42:57 2023

@author: sai
"""

class Node:
    def __init__(self, val:int, nidx:int, nid:str, left=None, right=None):
        self.val=val
        self.nidx=nidx
        self.nid=nid
        self.left=left
        self.right=right

    def __str__(self) -> None:
        return f'{self.nidx} {self.nid} {self.val.count()} : \
              {self.left.nidx if(self.left is not None) else "-"} \
              {self.right.nidx if (self.right is not None) else "-"}'

    def to_csv(self, sep='\t') -> None:
        """
        Returns the Node index, Node Id (Species ID), sketch count, load factor
        and indices of the left and right child nodes 

        Parameters
        ----------
        sep : str, optional
            DESCRIPTION. The default is '\t'.

        Returns
        -------
        None
            DESCRIPTION.

        """
        val_cnt=self.val.count()
        load_factor=round((val_cnt[0]/sum(val_cnt))*100.0,2)
        return f'{self.nidx}{sep}{self.nid}{sep}{self.val.count()}{sep}{load_factor}{sep}\
              {self.left.nidx if(self.left is not None) else "-"}{sep}\
              {self.right.nidx if (self.right is not None) else "-"}'
