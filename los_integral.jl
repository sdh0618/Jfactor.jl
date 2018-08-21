using Distributed

function nuclear_bulge(M)
	d2r = pi/180
	Rsun = 8.25
       	b=M[1]
	l=M[2]
	z = d*sin(b*d2r)
       	x = Rsun - d*cos(b*d2r)*cos(l*d2r)
       	y = -d*cos(b*d2r)*sin(l*d2r)
       	R = sqrt(Rsun^2-2*d*Rsun*cos(b*d2r)*cos(l*d2r)+d^2)
       	r = sqrt(Rsun^2-2*d*Rsun*cos(b*d2r)*cos(l*d2r)+(d*cos(b*d2r))^2)
       	if R < 0.006
       	nsc = 3.3e+15/(1+(R/2.2e-4)^2)
       	elseif 0.006 <= R <= 0.23
       	nsc = 27.2375*3.3e+15/(1+(R/2.2e-4)^3)
       	else
       	nsc = 0
       	end
        if r < 0.12
		nsd = 3e+11*(1000*r)^(-0.1)*exp(-abs(z)/0.045)
	elseif 0.12 <= r < 0.22
		nsd = 1.17278e+7*3e+11*(1000*r)^(-3.5)*exp(-abs(z)/0.045)
	elseif 0.22 <= r < 0.23
		nsd = 1.97226e+22*3e+11*(1000*r)^(-10)*exp(-abs(z)/0.045)
	else
		nsd = 0
	end

	return cos(b)*(nsc + nsd)
end

A = [ [b,l] for b=-89.75:0.5:89.75, l=-179.75:0.5:179.75 ]

jfactor =  @distributed (+) for d = 8:0.01:8.5
	pmap(nuclear_bulge,A)
end
