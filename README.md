# Compare different quadrature rules for integration

There are two examples provided for calculating the weights and abscissas for gaussian quadrature rules, try:

```
make
./gqconstants
```

or

```
python gqconstants.py
```

You can also use the C++ example as a guide to build your own executable

There is no need to look at rules >~25 for Gaussian quadrature.  And you can also stop at ~ 1000 divisions for the trapezoidal and Simpson's rules.  If you run much longer you'll see the numerical errors bevome visible for the trapezoidal, but hyou'll need to think about how to code efficiently or the running time may be very long.

