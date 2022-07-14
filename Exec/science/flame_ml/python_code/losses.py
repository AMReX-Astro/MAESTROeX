import torch
import torch.nn as nn

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def my_heaviside(x, input):
    y = torch.zeros_like(x)
    y[x < 0] = 0
    y[x > 0] = 1
    y[x == 0] = input
    return y

def my_exp_weights(x, b=6.0):
    w = torch.zeros_like(x, requires_grad=False)
    w[x < 0.25] = torch.pow(10., -b/(1 + torch.exp(-100*x[x < 0.25] + 5)) + b)
    w[x >= 0.25] = torch.pow(10., b/(1 + torch.exp(-100*x[x >= 0.25] + 45)))
    return w

def component_loss_f(prediction, targets):
    #Takes the MSE of each component and returns the array of losses
    loss = torch.zeros(prediction.shape[1])
    L = nn.MSELoss()
    for i in range(prediction.shape[1]):
        loss[i] = L(prediction[:, i], targets[:,i])
    return loss


def component_loss_f_L1(prediction, targets):
    #Takes the MSE of each component and returns the array of losses
    loss = torch.zeros(prediction.shape[1])
    L = nn.L1Loss()
    for i in range(prediction.shape[1]):
        loss[i] = L(prediction[:, i], targets[:,i])
    return loss


def logX_loss(prediction, target, nnuc=13):
    # We are working with mass fractions in the form of -1/log(X_k)

    X = prediction[:, :nnuc]
    X_target = target[:, :nnuc]
    enuc = prediction[:, nnuc]
    enuc_target = target[:, nnuc]

    L = nn.MSELoss()
    F = nn.L1Loss()

    # enuc is allowed to be negative
    # but penalty should be given if prediction is of different signs
    enuc_fac = 1
    enuc_loss = L(enuc, enuc_target) + enuc_fac * F(torch.sign(enuc), torch.sign(enuc_target))

    # we do not want negative values for mass fractions
    if torch.sum(X < 0) > 0:
        #how much do we hate negative numbers?
        factor = 1000  #a lot
    else:
        factor = 1

    return factor * L(X, X_target) + enuc_loss

def logX_loss_noenuc(prediction, target, nnuc=2):
    # We are working with mass fractions in the form of -1/log(X_k)

    X = prediction[:, :nnuc]
    X_target = target[:, :nnuc]

    L = nn.MSELoss()

    # we do not want negative values for mass fractions
    if torch.sum(X < 0) > 0:
        #how much do we hate negative numbers?
        factor = 1000  #a lot
    else:
        factor = 1

    return factor * L(X, X_target)

def loss_wexp_noenuc(prediction, target, nnuc=2, offset=6.0):
    X = prediction[:, :nnuc]
    X_target = target[:, :nnuc]

    L = nn.MSELoss()

    # use weights that are inversely proportional to the mass fractions
    wexp = my_exp_weights(X_target, offset)

    return torch.sum(wexp * L(X, X_target))

def rms_weighted_error(input, target, solution, atol=1e-6, rtol=1e-6):
    error_weight = atol + rtol * torch.abs(solution)
    weighted_error = (input - target) / error_weight
    rms_weighted_error = torch.sqrt((weighted_error**2).sum() / input.data.nelement())
    return rms_weighted_error


def loss_mass_fraction_conserv(prediction, target, nnuc=2):
    L = nn.MSELoss()
    total_pred = torch.sum(prediction[:, :nnuc], 1)
    total = torch.sum(target[:, :nnuc], 1)
    total.requires_grad = False

    return L(total_pred, total)

def loss_mass_fraction_log_conserv(factor, prediction, target, nnuc=2):
    L = nn.MSELoss()

    Xk = torch.exp(-factor/target[:, :nnuc])
    total = torch.sum(Xk, 1)
    total.requires_grad = False

    Xk_pred = torch.exp(-factor/prediction[:, :nnuc])
    total_pred = torch.sum(Xk_pred, 1)

    return L(total_pred, total)
