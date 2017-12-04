% Numerical solution of the 2D incompressible convective Magnetohydrodynamics 
% equations on a Square Domain [0,1]x[0,1] using a Fourier pseudo-spectral
% method and Implicit midpoint rule timestepping. 
%
%Periodic free-slip boundary conditions and Initial conditions:
    %u(x,z,0)=sin(2*pi*x)cos(2*pi*z)
    %v(x,z,0)=-cos(2*pi*x)sin(2*pi*z)
    %theta(x,y,0)=exp(-2.0*((4*pi*(xx-0.25)).^2+(4*pi*(zz-0.25)).^2));
% internal heating function
% sin(4*pi*xx+4*pi*zz)
%
% Draws on work on "Infinte prandtl convection on the 2D torus" by Doering, Muite and Whitehead (forthcoming)
clear all; format compact; format short; clc; clf;
 
Re=100;%Reynolds number
Rem=100; % Magnetic Reynolds number
Ra=1; % Rayleigh number
Dtheta=0.01;% Heat diffusion constant
%grid
Nx=64; h=1/Nx; x=h*(1:Nx);
Nz=64; h=1/Nz; z=h*(1:Nz)';
[xx,zz]=meshgrid(x,z);
 
%initial conditions for velocity field
u=sin(2*pi*xx).*cos(2*pi*zz);
v=-cos(2*pi*xx).*sin(2*pi*zz);
u_z=-2*pi*sin(2*pi*xx).*sin(2*pi*zz);
v_x=2*pi*sin(2*pi*xx).*sin(2*pi*zz);
omega=v_x-u_z;
% initial magnetic potential field
alpha=sin(4*pi*xx).*cos(6*pi*zz);
% initial temperature field
theta=exp(-2.0*((4*pi*(xx-0.25)).^2+(4*pi*(zz-0.25)).^2));
% internal heating function
intheat=sin(4*pi*xx+4*pi*zz);

% timestepping parameters
dt=0.0025; t(1)=0; tmax=1.0;
nplots=ceil(tmax/dt);
 
%wave numbers for derivatives
k_x=2*pi*(1i*[(0:Nx/2-1)  0 1-Nx/2:-1]');
k_z=2*pi*(1i*[(0:Nz/2-1)  0 1-Nz/2:-1]);
k2x=k_x.^2;
k2z=k_z.^2;
 
%wave number grid for multiplying matricies
[kxx,kzz]=meshgrid(k2x,k2z);
[kx,kz]=meshgrid(k_x,k_z);
 
% use a high tolerance so time stepping errors
% are not dominated by errors in solution to nonlinear
% system
tol=10^(-10);
%compute \hat{Phi}^{n+1,k+1}  
alphahat=fft2(alpha);
phihat=-alphahat./(kxx+kzz);  
 
%NOTE: kxx+kzz has to be zero at the following points to avoid a
% discontinuity. However, we suppose that the streamfunction has
% mean value zero, so we set them equal to zero
phihat(1,1)=0;                  
phihat(Nx/2+1,Nz/2+1)=0;
phihat(Nx/2+1,1)=0;
phihat(1,Nz/2+1)=0;
 
%computes {\psi}_x by differentiation via FFT
dphix = real(ifft2(phihat.*kx));  
%computes {\psi}_z by differentiation via FFT
dphiz = real(ifft2(phihat.*kz));
% components of magnetic field
bx=dphiz;    
bz=-dphix;   
 
%compute \hat{\omega}^{n+1,k}
omegahat=fft2(omega); alphahat=fft2(alpha);
intheathat=fft2(intheat);
thetatotal(1)=sum(sum(theta.^2))*h*h;
for i=1:nplots
    chg=1;
    % save old values
    uold=u; vold=v;  omegaold=omega; omegacheck=omega;
    omegahatold=omegahat; thetaold=theta; thetacheck=theta;
    alphahatold=alphahat; alphaold=alpha; alphacheck=alpha;
    bxold=bx; bzold=bz;
    while chg>tol
        % Fluid field
        % nonlinear {n+1,k}
        nonlinhat=0.25*fft2((u+uold).*ifft2((omegahat+omegahatold).*kx)+...
                            (v+vold).*ifft2((omegahat+omegahatold).*kz)-...
                            (bx+bxold).*ifft2((alphahat+alphahatold).*kx)-...
                            (bz+bzold).*ifft2((alphahat+alphahatold).*kz));
 
        %Implicit midpoint rule timestepping 
        omegahat=((1/dt + 0.5*(1/Re)*(kxx+kzz)).*omegahatold ...
                  +0.5*(Ra/Re)*kx.*(theta+thetaold)-nonlinhat)...
            ./(1/dt -0.5*(1/Re)*(kxx+kzz));
 
        %compute \hat{\psi}^{n+1,k+1}    
        psihat=-omegahat./(kxx+kzz);  
 
        %NOTE: kxx+kzz has to be zero at the following points to avoid a
        % discontinuity. However, we suppose that the streamfunction has
        % mean value zero, so we set them equal to zero
        psihat(1,1)=0;                  
        psihat(Nx/2+1,Nz/2+1)=0;
        psihat(Nx/2+1,1)=0;
        psihat(1,Nz/2+1)=0;
 
        %computes {\psi}_x by differentiation via FFT
        dpsix = real(ifft2(psihat.*kx));  
        %computes {\psi}_y by differentiation via FFT
        dpsiz = real(ifft2(psihat.*kz));
 
        u=dpsiz;    %u^{n+1,k+1}
        v=-dpsix;   %v^{n+1,k+1}
 
        % magnetic field
        nonlinhat=0.25*fft2((u+uold).*ifft2((alphahat+alphahatold).*kx)+...
                            (v+vold).*ifft2((alphahat+alphahatold).*kz)-...
                            (bx+bxold).*ifft2((omegahat+omegahatold).*kx)-...
                            (bz+bzold).*ifft2((omegahat+omegahatold).*kz));
 
        %Implicit midpoint rule timestepping 
        alphahat=((1/dt + 0.5*(1/Re)*(kxx+kzz)).*alphahatold-nonlinhat)...
            ./(1/dt -0.5*(1/Rem)*(kxx+kzz));
 
        %compute \hat{\psi}^{n+1,k+1}    
        phihat=-alphahat./(kxx+kzz);  
 
        %NOTE: kxx+kyy has to be zero at the following points to avoid a
        % discontinuity. However, we suppose that the streamfunction has
        % mean value zero, so we set them equal to zero
        phihat(1,1)=0;                  
        phihat(Nx/2+1,Nz/2+1)=0;
        phihat(Nx/2+1,1)=0;
        phihat(1,Nz/2+1)=0;
 
        %computes {\psi}_x by differentiation via FFT
        dphix = real(ifft2(phihat.*kx));  
        %computes {\psi}_z by differentiation via FFT
        dphiz = real(ifft2(phihat.*kz));
 
        bx=dphiz;    %u^{n+1,k+1}
        bz=-dphix;   %v^{n+1,k+1}
 
 
        thetax=0.5*ifft2(kx.*fft2(theta+thetaold));
        thetaz=0.5*ifft2(kz.*fft2(theta+thetaold));
        theta=ifft2((fft2(thetaold-dt*0.5*((uold+u).*thetax+(vold+v).*thetaz))+...
                    dt*0.5*Dtheta*(kxx+kzz).*fft2(thetaold) + dt*intheathat ...
                    )./(1-dt*0.5*Dtheta*(kxx+kzz)));
 
        %\omega^{n+1,k+1}
        omega=ifft2(omegahat);
        % check for convergence
        chg=max(max(abs(omega-omegacheck))) + max(max(abs(theta-thetacheck)))+...
            max(max(abs(alpha-alphacheck)))
 
        % store omega and theta to check for convergence of next iteration
        omegacheck=omega; thetacheck=theta; alphacheck=alpha;
    end     
     t(i+1)=t(i)+dt;  
     thetatotal(i+1)=sum(sum(theta.^2))*h*h;
     figure(1); 
     subplot(2,2,1);
     pcolor(xx,zz,real(omega));  shading interp; xlabel x; ylabel z; 
     title(['Fluid Vorticity, Time ',num2str(t(i+1))]); colorbar; drawnow;
     subplot(2,2,2);
     pcolor(xx,zz,alpha); shading interp;  xlabel x; ylabel z; 
     title(['Alpha, Time ',num2str(t(i+1))]); colorbar; drawnow;
     subplot(2,2,3);
     pcolor(xx,zz,real(ifft2(psihat))); shading interp; xlabel x; ylabel z; 
     title(['Psi, Time ',num2str(t(i+1))]); colorbar; drawnow;
     subplot(2,2,4);
     pcolor(xx,zz,real(theta)); shading interp;  xlabel x; ylabel z; 
     title(['Theta, Time ',num2str(t(i+1))]); colorbar; drawnow;
end
