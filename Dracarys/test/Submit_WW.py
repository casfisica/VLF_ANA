import commands as cmd
import time

ListOfSamples=['/WW_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM']

ITime=time.localtime()
StringITime=str(ITime.tm_year)+'_'+str(ITime.tm_mon)+'_'+str(ITime.tm_mday)+'_'+str(ITime.tm_hour)+'_'+str(ITime.tm_min)
print "Submitting Jobs "+StringITime

WorkingFolder='CrabConfigs_'+StringITime
cmd.getoutput('mkdir '+WorkingFolder)

for i in ListOfSamples:
    UsefulString=i.split('/')[1]+"_"+i.split('/')[2]
    cmd.getoutput('cp Crab_Template.py '+WorkingFolder+'/Crab_'+UsefulString+'.py')
    cmd.getoutput('sed -i -- "s/OUTPUTDIR/WW/g" '+WorkingFolder+'/Crab_'+UsefulString+'.py')
    cmd.getoutput('sed -i -- "s/TASK/'+UsefulString.split('Tune')[0]+i.split('/')[2].split("_")[-1]+'_'+StringITime+'/g" '+WorkingFolder+'/Crab_'+UsefulString+'.py')
    cmd.getoutput('sed -i -- "s#DATASAMPLE#'+i+'#g" '+WorkingFolder+'/Crab_'+UsefulString+'.py')
    print "Job being submitted: "+UsefulString
    CrabOutput=cmd.getoutput('crab submit -c '+WorkingFolder+'/Crab_'+UsefulString+'.py')
    print CrabOutput
    print "---------------------------------------------------------------------------------------"
