import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm


def negloglik(B=None, sigma_motor=None, sigma_int=None, sigma_pert=None, 
              bias=None, sigma_pred=None, sigma_p=None, sigma_v=None, 
              sigma_comb=None, s=None, c=None, model=None, num_trials=None, 
              rotation=None, fit=True,  vis_fb=None, eta_p=None, 
              beta_p_sat=None, x_hand=None):
    if model == "pea":
        x_stl, _ = pea(sigma_int, B, bias, sigma_motor, num_trials, vis_fb, rotation, 
                    fit, x_hand=x_hand)
    elif model == "piece":
        x_stl, _ = piece(sigma_pert, sigma_comb, bias, sigma_motor, num_trials, 
                      vis_fb, rotation, fit, x_hand=x_hand)
    elif model == "ssm":
        x_stl, _ = ssm_ege(B, sigma_motor, rotation)
    elif model == "premo":
        x_stl, _ = premo(B, sigma_v, sigma_p, sigma_pred, eta_p, bias, sigma_motor, 
                      num_trials, vis_fb, rotation, fit, x_hand=x_hand)
    elif model == "rem":
        x_stl, _ = rem(sigma_comb, s, c, bias, sigma_motor, num_trials, vis_fb, rotation, 
                    fit, x_hand=x_hand)
    
    if model == "piece":
        sigma0 = 0.5  # some "typical" value for sigma_pert
        lam = 0.1  # regularization factor
        
        # Adding bias term to single-trial adaptation measure (mu) 
        nll = -np.sum(-0.5 * np.log(2 * np.pi * sigma_motor**2) - ((x_hand - (x_stl + bias))**2 / (2 * sigma_motor**2)))
        # nll = nll + lam * (sigma_pert - sigma0)**2

    else:
        # Adding bias term to single-trial adaptation measure (mu) 
        nll = -np.sum(-0.5 * np.log(2 * np.pi * sigma_motor**2) - ((x_hand - (x_stl + bias))**2 / (2 * sigma_motor**2)))
    
    return nll


def pea(sigma_int, B, bias, sigma_motor, num_trials, vis_fb, rotation, fit=False, x_hand=None):
    '''
    Returns:
        x_stl : state estimate
        x_hand : motor output
    '''
    T = 0
    x_stl = np.zeros(num_trials)
    
    if fit == True:
        x_hand = x_hand
    else:
        x_hand = np.zeros(num_trials)
        x_hand[0] = np.random.normal(0, sigma_motor)
    
    for i in range(num_trials - 1):
        if vis_fb[i] == 0:
            x_v = 0
            sigma_v = 1e2  # high value that won't cause numerical problems
        else:
            x_v = x_hand[i] + rotation[i]
            sigma_v = 1.179 + 0.384 * np.abs(x_v)  # from Zhang et al 
        
        # Compute estimated hand position
        w_int = (1 / sigma_int**2) / (1 / sigma_int**2 + 1 / sigma_v**2)
        w_v = (1 / sigma_v**2) / (1 / sigma_int**2 + 1 / sigma_v**2)
        xhat_hand = w_v * x_v 
        
        # Update rule
        x_stl[i + 1] = B * (T - xhat_hand) 
        
        if fit == False:
            x_hand[i + 1] = bias + x_stl[i + 1] + np.random.normal(0, sigma_motor)
        
    return x_stl, x_hand


def premo(B, sigma_v, sigma_p, sigma_pred, eta_p, bias, sigma_motor, 
          num_trials, vis_fb, rotation, fit=False, x_hand=None):
    '''
    Model parameters
        B
        sigma_v
        sigma_p
        sigma_pred
        eta_p
    '''

    T = 0
    beta_p_sat = 5  # hard coding to reduce number of free params
    sigma_u = sigma_pred
    x_stl = np.zeros(num_trials)

    if fit == True:
        x_hand = x_hand
    else:
        x_hand = np.zeros(num_trials)
        x_hand[0] = np.random.normal(0, sigma_motor)

    for i in range(num_trials - 1):
        # Visual cue
        if vis_fb[i] == 0:
            x_v = 0
            sigma_v = 1e2
        else:
            x_v = x_hand[i] + rotation[i]
        J_v = 1 / sigma_v**2

        # Proprioceptive cue
        x_p = x_hand[i]
        
        # Precision values
        J_p = 1 / sigma_p**2
        J_u = 1 / sigma_u**2
  
        # Weights for each modality
        w_v = J_v / (J_v + J_u)
        w_p = J_p / (J_p + J_u)
    
        # beta_p is proprio shift due to crossmodal recal from vision;
        # there could be a problem when x_v == x_p
        if x_v > x_p:
            beta_p = np.min([np.abs(beta_p_sat), 
                             np.abs(eta_p * (w_v * x_v - w_p * x_p))])
        else:    
            beta_p = -np.min([np.abs(beta_p_sat), 
                              np.abs(eta_p * (w_v * x_v - w_p * x_p))])
        
        # Perceived hand position     
        x_prop_per = w_p * x_p + beta_p  

        # Update rule
        x_stl[i + 1] = B * (T - x_prop_per)
        
        if fit == False:
            x_hand[i + 1] = bias + x_stl[i + 1] + np.random.normal(0, sigma_motor)
    
    return x_stl, x_hand


def rem(sigma_comb, s, c, bias, sigma_motor, num_trials, vis_fb, rotation, 
        fit=False, x_hand=None):
    '''
    Model parameters
        B
        sigma_Comb
        s 
        c
    '''
    
    T = 0
    x_stl = np.zeros(num_trials)
    
    if fit == True:
        x_hand = x_hand
    else:
        x_hand = np.zeros(num_trials)
        x_hand[0] = np.random.normal(0, sigma_motor)

    for i in range(num_trials - 1):
        # Proprioceptive cue (does not play a role)
        x_p = x_hand[i]

        # Visual cue
        if vis_fb[i] == 0:
            x_v = 0
            p_rel = 0
        else:
            x_v = x_hand[i] + rotation[i]
            p_rel = s * (np.exp(-(x_v)**2 / (2 * sigma_comb**2))) / (np.exp(-(x_v)**2 / (2 * sigma_comb**2)) + c)

        # Update rule
        x_stl[i + 1] = (T - p_rel * x_v)
        
        if fit == False:
            x_hand[i + 1] = bias + x_stl[i + 1] + np.random.normal(0, sigma_motor)
        
    return x_stl, x_hand


def piece(sigma_pert, sigma_comb, bias, sigma_motor, num_trials, vis_fb, 
          rotation, fit=False, x_hand=None):
    '''
    Model parameters
        sigma_comb: combination of sigma_prop and sigma_pred
        sigma_motor
         
    Returns
        x_state : state estimate
        x_f : motor output
    '''
    
    # Function for computing Gaussian log-probabilities
    f = lambda x, mu, sigma: -0.5 * np.log(2 * np.pi * sigma**2) - 0.5 * (x - mu)**2 / sigma**2
    
    # Possible endpoint locations
    x_grid = np.arange(-15, 15, 0.1)  

    # For vectorized code
    x_fs = x_grid.reshape((len(x_grid), 1))  # possible finger endpoint locations (col vec)
    d_xvs = x_grid.reshape((1, len(x_grid)))  # possible rotation sizes (row vec)
    x_vs = x_grid  # possible locations of visual cues
    x_ps = x_grid  # possible locations of proprioceptive cues

    # Ideal observer     
    x_f_hat = np.zeros(num_trials)
    x_state = np.zeros(num_trials)
    K = np.zeros(num_trials)
    mu_pert = 0
    prior_pert = 0.5
    
    # If fitting model to data, use actual hand data
    if fit == True:
        x_f = x_hand
    else: 
        x_f = np.zeros(num_trials) 
        x_f[0] = np.random.normal(0, sigma_motor)
    
    # Loop through trials
    for i in range(num_trials - 1):
        # Proprioceptive cue contains motor prediction (no way to dissociate)
        xhat_p = x_f[i]

        if vis_fb[i] == 0:
            xhat_v = 0
            sigma_v = 1e2
        else:
            xhat_v = x_f[i] + rotation[i]
            sigma_v = 1.179 + 0.384 * np.abs(xhat_v)  # from Zhang et al 
        J_v = 1 / sigma_v**2
        J_pert = 1 / sigma_pert**2
        K[i] = J_v / (J_v + J_pert)

        # Compute no perturbation likelihood (working with log-probs for numerical accuracy)
        loglik_nopert = (
            f(xhat_p, x_fs, sigma_comb) + f(xhat_v, x_fs, sigma_v) 
            + f(x_fs, bias, sigma_motor) 
        )
        loglik_nopert = logsumexp(loglik_nopert.flatten(), b=0.1)
        likelihood_nopert = np.exp(loglik_nopert)

        # Compute perturbation likelihood
        loglik_pert = (
            f(xhat_p, x_fs, sigma_comb) + f(xhat_v, x_fs + d_xvs, sigma_v) 
            + f(d_xvs, mu_pert, sigma_pert) + f(x_fs, bias, sigma_motor)  
        )
        loglik_pert = logsumexp(loglik_pert.flatten(), b=0.01)
        
        # To account for no-vis fb trials
        if vis_fb[i] == 1:
            likelihood_pert = np.exp(loglik_pert)
        else:
            likelihood_pert = 0

        # Posterior over Causal node
        normalization_const = prior_pert * likelihood_pert + ((1 - prior_pert) * likelihood_nopert)
        post_pert = (prior_pert * likelihood_pert) / normalization_const  # posterior over Cause

        # Simulate trial-by-trial adaptation
        x_state[i + 1] = post_pert * K[i] * (rotation[i]) * -1
        if fit == False: 
            x_f[i + 1] = bias + x_state[i + 1] + np.random.normal(0, sigma_motor)
        
    return x_state, x_f


def calc_bic(ll, num_params, num_trials):
    bic = -2 * ll + num_params * np.log(num_trials)
    
    return bic


def calc_aic(ll, num_params):
    aic = -2 * ll + 2 * num_params

    return aic
