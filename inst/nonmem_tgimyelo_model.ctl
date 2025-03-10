$PROB Run 1 Tumor growth inhibition
$INPUT C   TIME   AMT   DV   EVID  
$DATA tgimyelo.dat IGNORE=C
$SUBROUTINES ADVAN13
$MODEL NCOMPARTMENTS=11

$PK
lV = THETA(1)
lCirc0 = THETA(2)
lW0 = THETA(3)

lka = THETA(4)
lkel = THETA(5)
ll0 = THETA(6)
ll1 = THETA(7)
lk1 = THETA(8)

lEmax = THETA(9)
lEC50 = THETA(10)
lkprol = THETA(11)
lktr = THETA(12)
lkcirc = THETA(13)

V = lV * EXP(ETA(1))
Circ0 = lCirc0 * EXP(ETA(2))
W0 = lW0  * EXP(ETA(3))

ka = lka * EXP(ETA(4))
kel = lkel * EXP(ETA(5))
l0 = ll0 * EXP(ETA(6))
l1 = ll1 * EXP(ETA(7))
k1 = lk1 * EXP(ETA(8))

Emax = lEmax * EXP(ETA(9))
EC50 = lEC50 * EXP(ETA(10))
kprol = lkprol * EXP(ETA(11))
ktr = lktr * EXP(ETA(12))
kcirc = lkcirc * EXP(ETA(13))

A_0(3) = W0
A_0(7) = Circ0
A_0(8) = Circ0
A_0(9) = Circ0
A_0(10) = Circ0
A_0(11) = Circ0

$DES
DADT(1) = -ka*A(1)
DADT(2) = ka*A(1) - kel*A(2)

CC = A(2)/V
w = A(3)+A(4)+A(5)+A(6)
eff = 2*CC/(80+CC)

DADT(3) = (2*l0*l1*A(3)**2)/((l1+2*l0*A(3))*w) - eff*A(3)
DADT(4) = eff*A(3)-k1*A(4)
DADT(5) = k1*(A(4)-A(5))
DADT(6) = k1*(A(5)-A(6))

DADT(7) = kprol*A(7) * (1-Emax*CC/(EC50+CC)) * (Circ0/A(11))**0.1 - ktr*A(7)
DADT(8) = ktr*(A(7)-A(8))
DADT(9) = ktr*(A(8)-A(9))
DADT(10) = ktr*(A(9)-A(10))
DADT(11) = ktr*A(10) - kcirc*A(11)

$ERROR
Y = CC

$THETA
2
7
0.01

5
2
0.1
0.2
0.5

1
50
2
2
2

$OMEGA
0.01
0.01
0.01

0.01
0.01
0.01
0.01
0.01

0.01
0.01
0.01
0.01
0.01

$SIGMA
0.1

$EST METHOD=1 MAXEVAL=9999 INTER PRINT=5
$TABLE ID TIME DV IPRE=CIPRED AMT