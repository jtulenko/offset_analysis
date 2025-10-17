from scipy import stats
import numpy as np
import pylab as pl
import matplotlib.patches as patches


#this script calculates various Chi2 stats for a dataset of ages and errors, and (if necessary)
#removes extreme values until the distribution of ages is close enough to 'expected'.
#This particular one takes a user-definted numpy array as inputs and plots the original and
#'pruned' datasets as camelplots as well as some stats in plot title.

###########################################################################################
#                                                                                         #
# STEP 1: Calculate expected values (mu) for various distributions that might fit dataset #
#                                                                                         #
###########################################################################################

#function to calculate the expected value (mu) for a normal distribution
def norm_expected(ages,errors):
    mu_a = []
    mu_b = []

    for i in range(0,len(ages)):
        mu_ia = ages[i] / (errors[i] ** 2)
        mu_ib = 1 / (errors[i] ** 2)

        mu_a.append(mu_ia)
        mu_b.append(mu_ib)

    mu = np.sum(mu_a) / np.sum(mu_b)

    #print(mu)

    return mu

        
#######################################################################################
#                                                                                     #
# STEP 2: Calculate Chi-Squared statistic (X2) for individual values and full dataset #
#                                                                                     #
#######################################################################################

#function to calculate the X2 stat for each age in the dataset and the summary X2 stat for the dataset
#for a normal distribution model
def norm_X2(ages,errors):
    X2_list = []

    norm_mu = norm_expected(ages, errors)

    for i in range(0, len(ages)):
        X2_i = ((ages[i] - norm_mu) / errors[i]) ** 2

        X2_list.append(X2_i)

    X2 = np.sum(X2_list)

    #print(X2)

    return X2_list, X2


################################################
#                                              #
# STEP 3: Calculate a p-value for each dataset #
#                                              #
################################################

#as far as I can tell, the function would apply to all modeled distributions so 
#I'll just combine into one function
def pX2(ages, errors):
    _, X2_norm = norm_X2(ages,errors)
    #_, X2_lnorm = lognorm_X2(ages,errors)
    #_, X2_uni = uniform_X2(ages,errors)

    dof = len(ages) - 1

    p_norm = 1 - stats.chi2.cdf(X2_norm, dof)
    #p_lnorm = 1 - stats.chi2.cdf(X2_lnorm, dof)
    #p_uni = 1 - stats.chi2.cdf(X2_uni, dof)

    #print(p_norm)#, p_lnorm, p_uni)

    return p_norm #, p_lnorm, p_uni


##########################################################################
#                                                                        #
# STEP 4: Prune function by identifying age/error pair with highest X2_i #
#                                                                        #
##########################################################################

#the pruning function that calculates individual X2 values for each age/error and prunes the worst one.
def X2_pruning(ages,errors):
    p_value = pX2(ages, errors)

    X2_i, X2 = norm_X2(ages,errors)

    if p_value < 0.01:
        X2_iMax = np.argmax(X2_i)

        pruned_ages = np.delete(ages, X2_iMax)
        pruned_errors = np.delete(errors, X2_iMax)

        return True, pruned_ages, pruned_errors
    else:
        return False, ages, errors
    

#####################################################################################################
#                                                                                                   #
# STEP 5: Final step! Now go thru dataset, test Chi2 stat, and prune extreme values until Chi2 stat #
# looks good or the dataset is too small.                                                           #
#                                                                                                   #
#####################################################################################################

#here is the final function that folks will use.
def X2_iteration(ages,errors):
    ages_init = ages
    errors_init = errors

    ages_final = ages
    errors_final = errors

    while True and len(ages_final) >= (len(ages_init) // 2):

        prune, ages_final, errors_final = X2_pruning(ages_final, errors_final)

        if not prune:
            break

    _, X2_init = norm_X2(ages_init, errors_init)
    _, X2_final = norm_X2(ages_final, errors_final)
    pX2_init = pX2(ages_init, errors_init)
    pX2_final = pX2(ages_final, errors_final)

    
    sum_init = 0
    sum_final = 0
    x_fill_init = np.arange((norm_expected(ages_init, errors_init) - np.std(ages_init)), (norm_expected(ages_init, errors_init) + np.std(ages_init)), 50)
    x_fill_final = np.arange((norm_expected(ages_final, errors_final) - np.std(ages_final)), (norm_expected(ages_final, errors_final) + np.std(ages_final)), 50)
    
    if len(ages_final) >= (len(ages_init) // 2):
        prune_len = len(ages_init) - len(ages_final)
        title = f'Chi2 Stat: {round(pX2_final, 3)} after pruning {prune_len} from {len(ages_init)} total samples. Age: {round(norm_expected(ages_final, errors_final) / 1000, 1)} +/- {round(np.std(ages_final) / 1000, 1)}'

    else: 
        prune_len = len(ages_init) - len(ages_final)
        title = f'Null rejected, Chi2 Stat: {round(pX2_final, 3)} too low {prune_len} of {len(ages_init)} removed. Age undetermined'

    pl.ion()
    pl.figure(1)
    pl.xlim(((np.min(ages_init))-(np.mean(errors_init)*7)),((np.max(ages_init))+(np.mean(errors_init)*7)))
    pl.fill_between(x_fill_init, 0, 1, color = 'grey', alpha = 0.2)
    pl.fill_between(x_fill_final, 0, 1, color =(137/255, 207/255, 240/255), alpha = 0.5)
    pl.title(title)
    for i in range(0, len(ages_init)):
        x_init = np.arange(((np.min(ages_init))-(np.mean(errors_init)*5)),((np.max(ages_init))+(np.mean(errors_init)*5)),10)
        y_init = (1/(errors_init[i]*np.sqrt(2*np.pi)))*(np.exp(-0.5*(((x_init-ages_init[i])/errors_init[i])**2)))
        sum_init = sum_init + y_init

        pl.plot(x_init, y_init, 'r', alpha = 0.5)

    for i in range(0, len(ages_final)):
        x_final = np.arange(((np.min(ages_final))-(np.mean(errors_final)*5)),((np.max(ages_final))+(np.mean(errors_final)*5)),10)
        y_final = (1/(errors_final[i]*np.sqrt(2*np.pi)))*(np.exp(-0.5*(((x_final-ages_final[i])/errors_final[i])**2)))
        sum_final = sum_final + y_final

        pl.plot(x_final, y_final, 'g')

    pl.ylim(0,(np.max(sum_final) + (0.1*np.max(sum_final))))
    pl.plot(x_init, sum_init, color='grey', alpha = 0.5)
    pl.plot(x_final, sum_final, 'k')


    #print(ages_init, errors_init, X2_init, pX2_init)
    #print(ages_final, errors_final, X2_final, pX2_final)

    return ages_init, errors_init, ages_final, errors_final



#here, you can generate your own dataset and run the X2_iteration function to do full analysis and
#plot a nice fig.
#ages1 = np.array([21000, 21100, 21200, 21300, 23000, 25000])
#errors1 = np.array([400, 400, 400, 400, 400, 4000])

#result = X2_iteration(ages1,errors1)