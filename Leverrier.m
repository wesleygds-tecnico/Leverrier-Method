function [autoval,pol] = Leverrier( A )
n = length(A); % ordem da matriz
autoval = zeros(n,1); % inicializar vari�vel que vai conter os autovalores de A
s = zeros(n,1); % inicializando as somas
pol = zeros(n+1,1); % inicializando os coeficientes do polinomio caracter�stico
B = A; % matriz auxiliar
for i = 1:n
    s(i) = trace(B);
    B = A*B;
end
pol(1) = 1.;
% resolve o Teorema de Newton 7.1 dos slides
% lembrando que o �ndice que o Matlab ou Octave inicia � sempre 1 (e n�o 0)
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