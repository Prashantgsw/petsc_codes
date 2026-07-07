import numpy as np
import matplotlib.pyplot as plt

filename = input("Enter data file name: ")
data = np.loadtxt(filename)

x    = data[:,0]
rho  = data[:,1]
rhou = data[:,2]
rhoE = data[:,3]

gamma = 1.4
u = rhou / rho
p = (gamma - 1.0) * (rhoE - 0.5*rho*u*u)

t = float(input("Enter physical time for this snapshot: "))

# --- Sod exact solution (Toro's method) ---
x_disc = 0.5
rho_L, u_L, p_L = 1.0, 0.0, 1.0
rho_R, u_R, p_R = 0.125, 0.0, 0.1

def f_K(p, rhoK, pK, cK):
    if p > pK:
        AK = 2.0/((gamma+1.0)*rhoK)
        BK = (gamma-1.0)/(gamma+1.0)*pK
        return (p-pK)*np.sqrt(AK/(p+BK))
    else:
        return (2.0*cK/(gamma-1.0)) * ((p/pK)**((gamma-1.0)/(2.0*gamma)) - 1.0)

def f_K_deriv(p, rhoK, pK, cK):
    if p > pK:
        AK = 2.0/((gamma+1.0)*rhoK)
        BK = (gamma-1.0)/(gamma+1.0)*pK
        return np.sqrt(AK/(p+BK)) * (1.0 - (p-pK)/(2.0*(p+BK)))
    else:
        return (1.0/(rhoK*cK)) * (p/pK)**(-(gamma+1.0)/(2.0*gamma))

def solve_star_region():
    c_L = np.sqrt(gamma*p_L/rho_L)
    c_R = np.sqrt(gamma*p_R/rho_R)
    p_old = max(1e-10, 0.5*(p_L+p_R) - 0.125*(u_R-u_L)*(rho_L+rho_R)*(c_L+c_R))
    for _ in range(50):
        f_total = f_K(p_old, rho_L, p_L, c_L) + f_K(p_old, rho_R, p_R, c_R) + (u_R-u_L)
        f_deriv_total = f_K_deriv(p_old, rho_L, p_L, c_L) + f_K_deriv(p_old, rho_R, p_R, c_R)
        p_star = p_old - f_total/f_deriv_total
        p_star = max(p_star, 1e-10)
        if abs(p_star - p_old) < 1e-10:
            break
        p_old = p_star
    u_star = 0.5*(u_L+u_R) + 0.5*(f_K(p_star, rho_R, p_R, c_R) - f_K(p_star, rho_L, p_L, c_L))
    return p_star, u_star

def sample_sod(S, p_star, u_star):
    c_L = np.sqrt(gamma*p_L/rho_L)
    c_R = np.sqrt(gamma*p_R/rho_R)
    if S <= u_star:
        if p_star > p_L:
            rho_star = rho_L * ((p_star/p_L + (gamma-1)/(gamma+1)) /
                                 ((gamma-1)/(gamma+1)*p_star/p_L + 1))
            S_wave = u_L - c_L*np.sqrt((gamma+1)/(2*gamma)*p_star/p_L + (gamma-1)/(2*gamma))
            if S < S_wave:
                return rho_L, u_L, p_L
            else:
                return rho_star, u_star, p_star
        else:
            c_star = c_L*(p_star/p_L)**((gamma-1)/(2*gamma))
            S_head = u_L - c_L
            S_tail = u_star - c_star
            if S < S_head:
                return rho_L, u_L, p_L
            elif S > S_tail:
                return rho_L*(p_star/p_L)**(1/gamma), u_star, p_star
            else:
                u_loc = 2.0/(gamma+1) * (c_L + (gamma-1)/2*u_L + S)
                c_loc = 2.0/(gamma+1) * (c_L + (gamma-1)/2*(u_L - S))
                rho_loc = rho_L*(c_loc/c_L)**(2/(gamma-1))
                p_loc = p_L*(c_loc/c_L)**(2*gamma/(gamma-1))
                return rho_loc, u_loc, p_loc
    else:
        if p_star > p_R:
            rho_star = rho_R * ((p_star/p_R + (gamma-1)/(gamma+1)) /
                                 ((gamma-1)/(gamma+1)*p_star/p_R + 1))
            S_wave = u_R + c_R*np.sqrt((gamma+1)/(2*gamma)*p_star/p_R + (gamma-1)/(2*gamma))
            if S > S_wave:
                return rho_R, u_R, p_R
            else:
                return rho_star, u_star, p_star
        else:
            c_star = c_R*(p_star/p_R)**((gamma-1)/(2*gamma))
            S_head = u_R + c_R
            S_tail = u_star + c_star
            if S > S_head:
                return rho_R, u_R, p_R
            elif S < S_tail:
                return rho_R*(p_star/p_R)**(1/gamma), u_star, p_star
            else:
                u_loc = 2.0/(gamma+1) * (-c_R + (gamma-1)/2*u_R + S)
                c_loc = 2.0/(gamma+1) * (c_R - (gamma-1)/2*(u_R - S))
                rho_loc = rho_R*(c_loc/c_R)**(2/(gamma-1))
                p_loc = p_R*(c_loc/c_R)**(2*gamma/(gamma-1))
                return rho_loc, u_loc, p_loc

p_star, u_star = solve_star_region()

rho_exact = np.zeros_like(x)
u_exact   = np.zeros_like(x)
p_exact   = np.zeros_like(x)

for i, xp in enumerate(x):
    if t > 0:
        S = (xp - x_disc)/t
        rho_exact[i], u_exact[i], p_exact[i] = sample_sod(S, p_star, u_star)
    else:
        if xp < x_disc:
            rho_exact[i], u_exact[i], p_exact[i] = rho_L, u_L, p_L
        else:
            rho_exact[i], u_exact[i], p_exact[i] = rho_R, u_R, p_R

title_suffix = f" (t={t:.4f}, N={len(x)})"

plt.figure()
plt.plot(x, rho, label='numerical')
plt.plot(x, rho_exact, '--', label='exact')
plt.xlabel("x"); plt.ylabel(r"$\rho$")
plt.title("Density" + title_suffix)
plt.legend(); plt.grid(True)
plt.savefig("rho_compare.png")

plt.figure()
plt.plot(x, u, label='numerical')
plt.plot(x, u_exact, '--', label='exact')
plt.xlabel("x"); plt.ylabel(r"$u$")
plt.title("Velocity" + title_suffix)
plt.legend(); plt.grid(True)
plt.savefig("u_compare.png")

plt.figure()
plt.plot(x, p, label='numerical')
plt.plot(x, p_exact, '--', label='exact')
plt.xlabel("x"); plt.ylabel(r"$p$")
plt.title("Pressure" + title_suffix)
plt.legend(); plt.grid(True)
plt.savefig("p_compare.png")

plt.show()