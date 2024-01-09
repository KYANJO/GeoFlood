      integer maxregions,mregions
      real*8 xlowregion,ylowregion,xhiregion,yhiregion
      real*8 tlowregion,thiregion
      integer minlevelregion,maxlevelregion
      
      parameter (maxregions=2000)

      dimension xlowregion(maxregions)
      dimension ylowregion(maxregions)
      dimension xhiregion(maxregions)
      dimension yhiregion(maxregions)
      dimension minlevelregion(maxregions)
      dimension maxlevelregion(maxregions)
      dimension tlowregion(maxregions)
      dimension thiregion(maxregions)

      common /regionparams/ xlowregion,ylowregion,xhiregion,yhiregion,tlowregion,thiregion,minlevelregion,maxlevelregion,mregions