using LinearAlgebra, PyPlot
# 1D
function Hamiltonian(N, t1, t2;period=1)
    n = 2*N   
    hopping1 = zeros(ComplexF64,n-1)
    for i = 1:n-1
        if i % 2 == 0
            hopping1[i] = t2
        else
            hopping1[i] = t1
        end
    end
    H = diagm(1 => hopping1, -1=> hopping1)
    if period == 1
        H[1,end] = t2
        H[end,1] = t2
    end
    return H
end

N = 100
H = Hamiltonian(N,1,2;period=1)
E,V=eigen(H)
#figure()
#plot(E,"k*")
#grid()
#plot(V[:,N],"k.")

x = 0:2N-1
phi = diagm(exp.(im*2*pi*x/(2N-1)))

# Px
S = V[:,1:N]' * phi * V[:,1:N] 
imag.(log.(det(S)))/2/pi

# LDOS
sumV = 0
for i =1:N
    V0=abs.(V[:,i]).^2
    sumV = sumV .+ V0
end

sumVi = zeros(N)
for i = 1:N
    sumVi[i] = sum(sumV[(i-1)*2+1:i*2])
end
figure(figsize=(2.5,2.5));plot(sumVi,"k.");
ylim([0.5,1.5])
xlabel("site")
ylabel("LDOS")
grid()
