function [mu_hat,P_hat] = KFA(mu,S,z,Q,R,H,phi)
mu_bar=phi*mu;
P_bar=phi*S*phi'+Q;
K=P_bar*H'*inv(H*P_bar*H'+R);
mu_hat=mu_bar+K*(z'-H*mu_bar);
P_hat=P_bar-K*H*P_bar;
end