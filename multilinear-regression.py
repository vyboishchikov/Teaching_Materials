# Python code for multiple lineal regression Y = a*X1 + b*X2 + c
# in linear or logarithmic scale
import sklearn.linear_model, math

X1 = [ 0.50, 0.50, 1.00, 1.00]
X2 = [ 0.50, 1.00, 0.50, 1.00]
Y = [ 0.0477, 0.135, 0.382, 1.08]

LogX1 = True # Set LogX1 = False for a linear X1 scale
LogX2 = True # Set LogX2 = False for a linear X2 scale
LogY = True  # Set LogY = False for a linear Y scale

TextX1 = "X1"
TextX2 = "X2"
TextY = "Y"

# Logarithmic scale
if LogX1:
	X1 = [math.log(r) for r in X1]
	TextX1 = "ln"+TextX1
if LogX2:
	X2 = [math.log(r) for r in X2]
	TextX2 = "ln"+TextX2
if LogY:
	Y = [math.log(r) for r in Y]
	TextY = "ln"+TextY


X = [p for p in zip(X1,X2)]

regr = sklearn.linear_model.LinearRegression().fit(X, Y)

print(TextY,'= ',str(regr.coef_[0])+"*"+TextX1+' + '+str(regr.coef_[1])+"*"+TextX2,'+', regr.intercept_)
R2 = regr.score(X, Y)
print("R2 =%7.4f"%R2)
