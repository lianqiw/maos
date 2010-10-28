function [y]=MUV(x,RL)
    y=RL.M*x-(RL.U*(RL.V'*x));
