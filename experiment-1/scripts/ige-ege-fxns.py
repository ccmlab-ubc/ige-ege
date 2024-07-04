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
