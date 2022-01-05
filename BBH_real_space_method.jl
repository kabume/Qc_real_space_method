using LinearAlgebra, PyPlot
#Bloch BC
function Hamiltonian(kx,ky)
    gamma_x = 0.5  # hopping inside one unit cell
    lambda_x = 1   # hopping between unit cells
    gamma_y = gamma_x
    lambda_y = lambda_x
    h = zeros(ComplexF64, 4, 4)
    h[1, 3] = gamma_x + lambda_x * exp(im*kx)
    h[2, 4] = gamma_x + lambda_x * exp(-im*kx)
    h[1, 4] = gamma_y + lambda_y * exp(im*ky)
    h[2, 3] = -gamma_y - lambda_y * exp(-im*ky)
    h[3, 1] = conj(h[1, 3])
    h[4, 2] = conj(h[2, 4])
    h[4, 1] = conj(h[1, 4])
    h[3, 2] = conj(h[2, 3]) 
    return h
end

kx = -pi:0.05:pi
ky = -pi:0.05:pi
E = zeros(ComplexF64, length(kx), length(ky), 4)
i = 1
for kx in kx
    j = 1
    for ky in ky
        h = Hamiltonian(kx,ky)
        E[i,j,:] = eigvals(h)
        j += 1
    end
    i += 1
end
E = real.(E)
figure()
surf(E[:,:,1]);surf(E[:,:,2]);surf(E[:,:,3]);surf(E[:,:,4])

#OBC and CBC
function Hamiltonian(γ, λ, δ, N;period=1)
    n = 2N
    hopping1 = zeros(ComplexF64,n^2-1)
    for i = 1:n^2-1
        if i % n == 0
            hopping1[i] = 0
        elseif i % n % 2 == 0
            hopping1[i] = λ
        else
            hopping1[i] = γ
        end
    end

    hopping2 = zeros(ComplexF64,n^2-n)
    h2_1 = zeros(ComplexF64, n)
    for i = 1:n
        if i % 2 == 0
            h2_1[i] = γ
        else
            h2_1[i] = -γ 
        end
    end
    h2_2 = zeros(ComplexF64, n)
    for i = 1:n
        if i % 2 == 0
            h2_2[i] = λ
        else
            h2_2[i] = -λ 
        end
    end
    for i = 1:n-1
        if i % 2 == 0
            hopping2[n*(i-1)+1:n*(i-1)+n] = h2_2
        else
            hopping2[n*(i-1)+1:n*(i-1)+n] = h2_1
        end
    end
    onsite = zeros(ComplexF64, n^2)
    for i = 1:n^2
        if i % 2 == 0
            onsite[i] = -δ
        else
            onsite[i] = δ
        end
    end

    H = diagm(0 => onsite, 1 => hopping1, -1=> hopping1, n => hopping2, -n => hopping2)
    if period == 1 #PBC
        i1 = 1:n
        i2 = n^2-n+1:n^2
        for i in zip(i1,i2)
            if i[1] % 2 == 0
                H[i[1],i[2]] = -λ
                H[i[2],i[1]] = -λ
            else
                H[i[1],i[2]] = λ
                H[i[2],i[1]] = λ
            end
        end
        j1 = 1:n:n^2-n+1
        j2 = n:n:n^2
        for j in zip(j1,j2)
            H[j[1], j[2]] = λ
            H[j[2], j[1]] = λ
        end
    end
    return H
end

γ=0.01; λ = 1; δ = 0.1; N = 10
E,V = eigen(Hamiltonian(γ,λ,δ, N,period=0))
figure();plot(E,".")
figure();imshow(reshape(real.(V[:,2N^2+5]),2N,2N))

#LDOS
sumV = 0
for i =1:2*N^2
    norm = V[:,i]' * V[:,i]
    V0=abs.(V[:,i]).^2 / norm
    sumV = sumV .+ V0
end
sumV = reshape(sumV,2*N,2*N)
sumVi = zeros(N,N)
for i = 1:N
    for j =1:N
        sumVi[i,j] = sum(sumV[(i-1)*2+1:i*2, (j-1)*2+1:j*2])
    end
end
figure(figsize=(2.5,2.5));imshow(sumVi,cmap="coolwarm");colorbar()
clim([1.5,2.5])

#Qc: 0.5 or 0
x = 0:2N-1
y = copy(x)
phi = exp.(im*2*pi*x*y'/(2N-1)^2)
phi = diagm(phi[:])
S = V[:,1:2N^2]' * phi * V[:,1:2N^2] 
angle(det(S))/2/pi

#Px: keep zero because of C2 symmetry
x = 0:2N-1
y = ones(length(x))
phi = exp.(im*2pi*x*y'/(2N-1))
phi = diagm(phi[:])
S = V[:,1:2N^2]' * phi * V[:,1:2N^2] 
angle(det(S))/2/pi
