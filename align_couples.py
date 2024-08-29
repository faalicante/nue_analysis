import ROOT
import fedrarootlogon

brick=121
xmin=42400 - 2000
xmax=42400 + 2000
ymin=13000 - 2000
ymax=13000 + 2000

bin_size = 50
xbin = int((xmax-xmin)/bin_size)
ybin = int((ymax-ymin)/bin_size)

path='/eos/experiment/sndlhc/emulsionData/2022/emureco_Napoli/RUN1/b000121/sh_1/b000121'
histofile = ROOT.TFile(path+"/hist_couples_aligned.root","RECREATE")
ntuplefile = ROOT.TFile(path+"/ntuples_couples_aligned.root","RECREATE")
ntupletree = ROOT.TNtuple("segments","Ntuple of segments","p:x:y:tx:ty")
#scanset
sspath = path+"/.."
sproc = ROOT.EdbScanProc()
sproc.eProcDirClient=sspath
id = ROOT.EdbID(brick,0,0,0)
ss = sproc.ReadScanSet(id)
ss.Brick().SetID(brick)

npl = ss.eIDS.GetEntries()
cutstring = f"eCHI2P<2.4&&s.eW>20&&eN1<=1&&eN2<=1&&s1.eFlag>=0&&s2.eFlag>=0"
cut = ROOT.TCut(cutstring)
cut.Print()
for i in range(npl):
  idplate = ss.GetID(i)
  nplate = idplate.ePlate
  plate = ss.GetPlate(idplate.ePlate)
  p = ROOT.EdbPattern()
  ect = ROOT.EdbCouplesTree()
  ect.InitCouplesTree("couples",path+f"/p{nplate:03}/{brick}.{nplate}.0.0.cp.root","READ")
  #addingcut
  ect.eCut = cut 
  cutlist = ect.InitCutList()
  if not ect.eTree: continue #no TTree in file
  if not cutlist: 
    print("We have no entries, quitting!")
    continue
  nlst =cutlist.GetN()
  print("We have", nlst, "good couples in plate", nplate)
  
  nsegcut = cutlist.GetN()
  nseg = ect.eTree.GetEntries()
  hxy = ROOT.TH2D(f"h_{nplate}",f"couples_plate_{nplate};x[#mum];y[#mum]", xbin, xmin, xmax, ybin, ymin, ymax)
  for ientry in range(nsegcut):
    iseg = cutlist.GetEntry(ientry)
    ect.GetEntry(iseg)
    seg=ect.eS
    #setting z and affine transformation
    seg.SetZ(plate.Z())
    seg.SetPID(i)
    seg.Transform(plate.GetAffineXY())
    seg.Transform(plate.GetAffineTXTY())
    hxy.Fill(seg.X(), seg.Y())
    ntupletree.Fill(nplate, seg.X(), seg.Y(), seg.TX(), seg.TY())
  histofile.cd()
  hxy.Write()
ect.Close()
histofile.Close()



