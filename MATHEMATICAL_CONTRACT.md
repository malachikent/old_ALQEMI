# MATHEMATICAL CONTRACT:

$$

[\hat{A},\hat{B_{1}}\hat{B_{2}}] = [\hat{A},\hat{B_{1}}]\hat{B_{2}} + \hat{B_{1}}[\hat{A},\hat{B_{2}}]
$$

$$

[\hat{A},\hat{B_{1}} … \hat{B_{n}}] = \sum_{k=1}^{n} \hat{B_{1}} … \hat{B_{k-1}}[\hat{A},\hat{B_{k}}]\hat{B_{k+1}}…\hat{B_{n}}

$$


$$

[\hat{A},\hat{B_{1}}\hat{B_{2}}] = [\hat{A},\hat{B_{1}}]_+ \hat{B_{2}} - \hat{B_{1}}[\hat{A},\hat{B_{2}}]_+

$$

$$

[\hat{A},\hat{B_{1}} … \hat{B_{n}}] = \sum_{k=1}^{n} (-1)^{k-1} \hat{B_{1}} … [\hat{A},\hat{B_{k}}]_+ … \hat{B_{n}}

$$ (for n= even number)

$$

[\hat{A},\hat{B_{1}}\hat{B_{2}}]_+ = [\hat{A},\hat{B_{1}}]\hat{B_{2}} + \hat{B_{1}}[\hat{A},\hat{B_{2}}]_+ =  [\hat{A},\hat{B_{1}}]_+ \hat{B_{2}} - \hat{B_{1}}[\hat{A},\hat{B_{2}}]

$$

$$

[\hat{A},\hat{B_{1}} … \hat{B_{n}}]_+ = \sum_{k=1}^{n} (-1)^{k-1} \hat{B_{1}} … [\hat{A},\hat{B_{k}}]_+ … \hat{B_{n}}

$$ (for n= odd number)

Examples:

$$

[a^{\dagger}_{P}, a^{\dagger}_{Q}a_{R}] = -\delta_{PR}a^{\dagger}_{Q}

$$

$$

[a_{P}, a^{\dagger}_{Q}a_{R}] = \delta_{PQ}a_{R}

$$

$$

[a^{\dagger}_{P}a_{Q}, a^{\dagger}_{R}a_{S}] = \delta_{QR}a^{\dagger}_{P}a_{S} - \delta_{PS}a^{\dagger}_{R}a_{Q}

$$

$$

[a^{\dagger}_{P}a_{Q}, a^{\dagger}_{R}a_{S}a^{\dagger}_{M}a_{N}] = \delta_{QR}a^{\dagger}_{P}a_{S}a^{\dagger}_{M}a_{N} - \delta_{PS}a^{\dagger}_{R}a_{Q}a^{\dagger}_{M}a_{N} + \delta_{QM}a^{\dagger}_{R}a_{S}a^{\dagger}_{P}a_{N} - \delta_{PN}a^{\dagger}_{R}a_{S}a^{\dagger}_{M}a_{Q}
$$

$$
[a^{\dagger}_{P}a_{Q}, [a^{\dagger}_{R}a_{S},a^{\dagger}_{M}a_{N}]] = \delta_{SM}\delta_{QR}a^{\dagger}_{P}a_{N} - \delta_{SM}\delta_{PN}a^{\dagger}_{R}a_{Q} - \delta_{RN}\delta_{QM}a^{\dagger}_{P}a_{S} + \delta_{RN}\delta_{PS}a^{\dagger}_{M}a_{Q}
$$
