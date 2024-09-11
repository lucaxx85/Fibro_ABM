classdef ImmuneStates < uint32
    enumeration
        M0Static      (1)
        Empty         (2)
        M0Moving      (3)
        M1Static      (4)
        M1Moving      (5)
        M2Static      (6)
        M2Moving      (7)
        MIntStatic    (8) % intermediate: between M1/M2
        MIntMoving    (9) % intermediate: between M1/M2
        F0Static       (10)
        F0Moving       (11)
        F1Static       (12)
        F1Moving       (13)
    end
end