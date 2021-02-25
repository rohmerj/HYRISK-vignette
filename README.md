# HYRISK-vignette
This is the vignette of the R package [HYRISK](www). Click [here](https://rawcdn.githack.com/rohmerj/HYRISK-vignette/aaa78cc794e8cc0ff694f7d89356b7b532219d3a/hyrisk_demo.html)

```{r}
library(HYRISK)
```

## Description of the case

The case study is focused on the stability analysis of a dyke described by [@Ferson06]. The dyke has revetments made of masonry blocks subject to wave action as depicted schematically in the Figure below. The stability is estimated as the difference between the dike strength minus the stress acting on it as $Z = strength - stress =\Delta.D -\frac{H.\tan(\alpha)}{\cos(\alpha).M.s^{0.5}}$ where $\Delta$ is the relative density of the revetment blocks, $D$ is their thickness, $\alpha$ is the slope of the revetment. The wave characteristics are the significant wave height $H$, and the offshore peak wave steepness $s$. The factor $M$ reflects the risk analyst$^{'}$s vision on the uncertainty related to the model itself, i.e. its ability to reproduce reality. 

If $Z\ge0$, the dike is stable (the strength is greater than the stress); unstable otherwise. 
The study is focused on the estimate of the probability for $Z$ to become negative, which is considered a measure of the dike reliability.

![A schematic overview of the dyke.](images/dyke.png)

The assessment function is coded as follows:
```{r}
FUN<-function(X){
  delta=X[1]
  D=X[2]
  alpha=X[3]
  M=X[4]
  H=X[5]
  s=X[6]
return(delta*D-(H*tan(alpha)/(cos(alpha)*M*sqrt(s))))
}

```
## Step 1: uncertainty representation

The first step focuses on uncertainty representation. 
It aims at selecting the most appropriate mathematical tool to represent uncertainty on the considered parameter. 
The available options are: interval, possibility distribution (trapezoidal or triangular, see e.g. \cite{Baudrit06a}), probability distribution (normal, lognormal, beta, triangle, uniform, Gumbel or user-defined), 
a probability distribution with imprecise parameters, i.e. a family of parametric probability distributions represented by a p-box [@Ferson02]. For the sake of clarity, we use the generic term imprecise probability 
to designate such a uncertainty representation tool. 
The procedure in *HYRISK* first uses the *CREATE_INPUT* function to define the input variables (imprecise, random or fixed); for instance by setting the values of the bounds of the interval, the mean and the standard deviation of a normal probability distribution, etc. Second, the *CREATE_DISTR* function assigns the corresponding distribution (probability or possibility) to each uncertain input. 


Parameter          | Symbol | Uncertainty type | Representation
-------------------|--------|------------------|-----------------
Significant wave height| $H$ | Randomness | p-box of type Weibull with imprecise shape and scale
Weibull scale         | $\lambda$ | Imprecision | Interval $[1.2, 1.5]$ 
Weibull shape          | $k$ | Imprecision | Interval $[10, 12 ]$ 
Peak wave steepness          | $s$ | Randomness | Normal (Gaussian) probability distribution with mean=0.040 and standard deviation=0.0055 
Revetment density          | $\Delta$ | Imprecision | Interval $[1.60, 1.65]$ 
Revetment thickness         | $H$ | Imprecision | Interval $[0.68, 0.72]$ 
Slope angle          | $\lambda$ | Imprecision | Triangular possibility distribution of support $[0.309, 0.328]$ and core ${0.318}$
Expert-defined factor         | $M$ | Imprecision | Trapezoidal possibility distribution of support $[3, 5.2]$ and core $[4, 5]$. 

```{r}
#INPUT PARAMETERS
ninput<-8 
#Number of input parameters 
# input parameters for the model 
# + subvariables in case of a 2-level propagation
input<-vector(mode="list", length=ninput)

input[[1]]=CREATE_INPUT(
                        name=expression(Delta),
                        type="possi",distr="interval",param=c(1.6,1.65))
input[[2]]=CREATE_INPUT(
                        name="D",
                        type="possi",distr="interval",param=c(0.68,0.72))
input[[3]]=CREATE_INPUT(
                        name=expression(alpha),
                        type="possi",distr="triangle",
                        param=c(0.309,0.318,0.328))
input[[4]]=CREATE_INPUT(
                        name="M",
                        type="possi",distr="trapeze",param=c(3,4,5,5.2))
input[[5]]=CREATE_INPUT(
                        name="H",
                        type="impr proba",distr="user",param=c(7,8),
                        quser=qweibull,ruser=rweibull)
input[[6]]=CREATE_INPUT(
                        name="s",
                        type="proba",distr="normal",param=c(0.03,0.006))
input[[7]]=CREATE_INPUT(
                        name="k",
                        type="possi",distr="interval",param=c(10,12))
input[[8]]=CREATE_INPUT(
                        name=expression(lambda),
                        type="possi",distr="interval",param=c(1.2,1.5))

#CREATION OF THE DISTRIBUTIONS ASSOCIATED TO THE PARAMETERS
input=CREATE_DISTR(input)
```
A visualisation function *PLOT_INPUT* has also been implemented to handle the different representation forms. 
```{r}
#VISUALISATION of INPUT UNCERTAINTY REPRESENTATION
PLOT_INPUT(input)
```
![Uncertainty Representation for the dyke failure case.](images/representation.png)



