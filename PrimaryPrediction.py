import numpy as np
import matplotlib.pyplot as plt

#current state of hive
#Hstate={'hivebees':1000,'foragers':0,'eggs':0,'honey':0}
#Current global state
#Gstate={'time':0,'pollinated':0,'totalarea':445}

#function of proportion of max amount of eggs laid in a day given by t by the queen
def eggprop(t,ps=8):
    el=abs(np.sin(np.pi*t/365)**ps)
    return el

"""
colony_sim simulates one colony's interactions in a day
E_hist - list of how many eggs in each day
H_hist - list of how many hivebees in each day
h,k - production of honey/nectar per acre (g/acre)
R - rate of acre pollination per day (acre/bee*day)
"""
def colony_sim(E_hist,H_hist,h,k,R,hive_state,global_state,Dh,o,df,hs,d,b,v,L,Hm,w):
    E_past=(E_hist[-hs-1]if len(E_hist)>=hs+1 else 0)
    H_past=(H_hist[-d-1]if len(H_hist)>=d+1 else 0)
    dP=R*hive_state['foragers']*(global_state['totalarea']-global_state['pollinated'])/global_state['totalarea']
    Hs=hive_state['honey']+h*dP*((Hm-hive_state['honey'])/Hm)-(E_past*b*(1-np.exp(-hive_state['honey']/(E_past*b))) if E_past>0 else 0)
    f=Hs+k*dP
    B=hive_state['hivebees']+hive_state['foragers']
    S=(o*(hive_state['foragers']**2/B) if B>0 else 0)
    dHB=(E_past*(f**2/(f**2+(b*E_past)**2))*(hive_state['hivebees']/(hive_state['hivebees']+v)) if E_past!=0 else 0)+S
    dHD=-H_past
    dFB=H_past*Dh
    dFD=-S-df*(hive_state['foragers']-S)
    dEB=L*eggprop(global_state['time'])*(B/(B+w))
    dED=-E_past
    H=hive_state['hivebees']+dHB+dHD
    F=hive_state['foragers']+dFB+dFD
    E=hive_state['eggs']+dEB+dED
    H_hist.append(dHB)
    E_hist.append(dEB)
    P=global_state['pollinated']+dP
    return {'pollinated':P,'hivebees':H,'foragers':F,'eggs':E,'honey':Hs}

'''
season function takes in parameters and simulates a few seasons, returning other important variables

IMPORTANT VARIABLES

t = first day of simulation
honeyprod - production of honey per acre (g/accre)
nectarprod - production of nectar per acre (g/acre)
rateacreprod - rate of acre pollination per day per bee (acre/bee*day)
MaxArea - max pollinatable area in acres
hivebees_in - initial honey bees
foragers_in - initial forager bees
eggs_in - initial eggs
honey_in - initial stored honey (g)
survive_hive - the proportion of hive bees that survive until becoming foragers
death_forag - the proportion of forager bees that die daily
seasons- how many seasons to simulate (default 2, but must be >2)

OUTPUTS

total pollinated land (acres)
boolean whether colony collapses
boolean whether colony population has decreased
boolean whether colony is stable (reaches equillibrium, whether it be with decreased pop)
boolean whether colony is still growing
'''

def season(t,honeyprod,nectarprod,rateacreprod, MaxArea, hivebees_in, foragers_in, eggs_in, honey_in,survive_hive=1,seasons=2,forag_revert=0.8,death_forag=0.0347,hatchtime=20,hivebeetime=4,necpollhalf=0.25,beehalf=2.5,maxegglay=3000,honeymax=5000,beepoplay=40000):
    if(seasons<2):
        raise ValueError("season parameter is less than 2")
    #memory from previous days
    timelist=[]
    beelist=[]
    #stores all initial variables into Hstate and Gstate which will be changed later on
    Hstate={'hivebees':hivebees_in,'foragers':foragers_in,'eggs':eggs_in,'honey':honey_in}
    Gstate={'time':t,'pollinated':0,'totalarea':MaxArea}
    seasont=0
    #memory from previous days
    Elist=[eggs_in]
    Hlist=[hivebees_in]
    totpoloutput=0
    while(seasont<365*seasons):
        Bcurrent=Hstate['hivebees']+Hstate['foragers']
        timelist.append(Gstate['time'])
        beelist.append(Bcurrent)
        daynew=colony_sim(Elist,Hlist,honeyprod,nectarprod,rateacreprod,Hstate,Gstate,survive_hive,forag_revert,death_forag,hatchtime,hivebeetime,necpollhalf,beehalf,maxegglay,honeymax,beepoplay)
        Hstate={'hivebees':daynew['hivebees'],'foragers':daynew['foragers'],'eggs':daynew['eggs'],'honey':daynew['honey']}
        Gstate={'time':Gstate['time']+1,'pollinated':(daynew['pollinated'] if (seasont+1)%365!=0 else 0),'totalarea':Gstate['totalarea']}
        if((seasont+1)%365==0):
            totpoloutput+=daynew['pollinated']
        seasont+=1
    initialmax=max(beelist[:364])
    penmax=max(beelist[-730:-366])
    finalmax=max(beelist[-365:])
    collapse=True if finalmax<0.01*initialmax else False
    diminished=True if finalmax<0.9*initialmax and not collapse else False
    stable=True if abs(finalmax-penmax)<0.1*penmax and not collapse else False
    increasing=True if finalmax>1.1*penmax else False
    return totpoloutput, collapse, diminished, stable, increasing

xlist=[]
ylist=[]
for i in range(15,65):
    ylist.append(season(0,1500,100,0.00002,445,170,120,0,0,necpollhalf=0.01*i,beehalf=0.1*i)[0])
    xlist.append(0.01*i-0.15)
plt.plot(xlist,ylist)
plt.xlabel("Death Rate")
plt.ylabel("Accumulated Pollinated Area over 2 years")
plt.title("Offset Egg Death Rate")
plt.show()