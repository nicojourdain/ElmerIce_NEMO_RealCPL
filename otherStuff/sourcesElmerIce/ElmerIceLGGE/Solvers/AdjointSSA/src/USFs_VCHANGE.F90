!#######################################################################
!#
!#  A collection of USER FUNCTIONS to perform variable changes in
!inverse methods; i.e. beta=10^a or beta=a^2
!#
!#######################################################################
!# Compute VarOut=10^VarIn
       FUNCTION TenPowerA(Model,nodenumber,VarIn) RESULT(VarOut)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn,VarOut

       VarOut = 10._dp**(VarIn)

       End FUNCTION TenPowerA
!# Compute DJDA from DJDB if B=10^A: DJDA=DJDB*ln(10)*10^A
!# DJDB=VarIn(1)
!# A=VarIn(2)
       FUNCTION Derivative_TenPowerA(Model,nodenumber,VarIn) RESULT(VarOut)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(2),VarOut

       VarOut = VarIn(1)*(10.0**(VarIn(2)))*log(10.0)

       End FUNCTION Derivative_TenPowerA
!# Compute VarOut=VarIn*VarIn
       FUNCTION Asquare(Model,nodenumber,VarIn) RESULT(VarOut)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn,VarOut

       VarOut = VarIn*VarIn
       END FUNCTION Asquare
!# Compute DJDA from DJDB if B=A^2: DJDA=DJDB*2A
!# DJDB=VarIn(1)
!# A=VarIn(2)
       FUNCTION Derivative_Asquare(Model,nodenumber,VarIn) RESULT(VarOut)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(2),VarOut

       VarOut = 2.0*VarIn(1)*VarIn(2)

       End FUNCTION Derivative_Asquare
!#
