import sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from certify_1d import certify
vals = {0.25:(0.6092595657,0.6029620602),0.5:(0.4556645645,0.63841333210),0.75:(0.3478935313,0.6706266),
        1.0:(0.2707000923,0.69988463),0.1:(0.7330166680,0.5800041202),0.0:(0.8330785997,0.5639490804)}
res = []
for a in sys.argv[1:]:
    k = float(a)
    zg = [-0.3,-0.25,-0.2,-0.15,-0.1] if k >= 1.0 else None
    res.append(certify(k, *vals[k], dalpha=1e-13, dq=1e-9, zgrid=zg))
json.dump(res, open(f"tight_{'_'.join(sys.argv[1:])}.json","w"), indent=1)
