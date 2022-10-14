function [autoval,pol] = Leverrier( A )
n = length(A); % ordem da matriz
autoval = zeros(n,1); % inicializar variável que vai conter os autovalores de A
s = zeros(n,1); % inicializando as somas
pol = zeros(n+1,1); % inicializando os coeficientes do polinomio característico
B = A; % matriz auxiliar
for i = 1:n
    s(i) = trace(B);
    B = A*B;
end
pol(1) = 1.;
% resolve o Teorema de Newton 7.1 dos slides
% lembrando que o índice que o Matlab ou Octave inicia é sempre 1 (e não 0)
for k = 2:n+1
    soma = 0.;
    for j = 1:k-1
        soma = soma + pol(j)*s(k-j);
    end
    soma;
    pol(k) = -soma/(k-1)
end
autoval = roots(pol);
end