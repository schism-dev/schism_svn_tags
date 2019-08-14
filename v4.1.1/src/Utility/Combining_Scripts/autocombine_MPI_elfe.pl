#! /usr/bin/perl -w

# This script can be run during or after a parallel SELFE run (e.g., using mpiexec or qsub);
# launch it shortly after the SELFE run starts to make sure outputs/local_to_global_* are there (for
# the same reason, purge old outputs/* before launching this script to get right # of CPUs etc).
# It is useful when you don't know the strucutre of the cluster.
# It runs as a daemon and will combine outputs along the way.
# Run on the run dir (i.e., one level above outputs/) on any system (but
# make sure the combine code below is compatible).
# Assume that *_elev.61 is always output, and that 
# the code will output (end_stack+1)_00* (otherwise the script will hang).
# For side-/center-based outputs, the combine code needs sidecenters.gr3/centers.gr3.

#use Cwd;

if(@ARGV != 2) 
{ 
  print "$0 <start # of stacks> <end # of stacks>\n";
  exit(1);
}
print "$0 @ARGV\n";
$start_stack=$ARGV[0]; $end_stack=$ARGV[1];  #$ntr=$ARGV[2]; $nproc=$ARGV[3];
$netcdf = (@ARGV == 3 ? $ARGV[2] : 0);

#Get nproc and ntracers
if(!-e "outputs") {die "No outputs dir!";}
open(FL,"<outputs/local_to_global_0000"); @all=<FL>; close(FL);
($_,$_,$_,$nproc,$ntr)=split(" ",$all[0]);

$code="~/bin/combine_output6";
@vars=('elev.61','hvel.64','pres.61','airt.61','shum.61','srad.61',
       'flsu.61','fllu.61','radu.61','radd.61','flux.61','evap.61',
       'prcp.61','bdrc.61','wind.62','wist.62','dahv.62','vert.63','temp.63',
       'salt.63','conc.63','tdff.63','vdff.63','kine.63','mixl.63',
       'zcor.63','qnon.63','totN.63','depth.61','qtot.62','qsus.62','qbdl.62',
       'dpdxy.62','bedd50.61','bstress.61','brough.61','cdsed.61',
       'cflsed.61','d50moy.61','qav',
       'hvel.67','vert.69','temp.70','salt.70','dtbe.66','z0cr.66',
       'z0sw.66','z0wr.66','z0st.66','z0eq.66','bpgr.65','wafo.67',
       'bdoc.66','bnh4.66','bno3.66','bpo4.66','bcod.66','sbdo.66','sbsa.66',
       'ptot.63');

for($i=1;$i<=$ntr;$i++)
{
  $name="bfrac_$i.61";
  @vars=(@vars,$name);
  $name="qbdl_$i.62";
  @vars=(@vars,$name);
  $name="trcr_$i.63";
  @vars=(@vars,$name);
} #for

for($i=1;$i<=25;$i++)
{
  if($i<=23)
  {$name="wwm_$i.61";}
  else
  {$name="wwm_$i.62";}
  @vars=(@vars,$name);
} #for

for($i=1;$i<=$ntr/2;$i++)
{
  $name="age_$i.63";
  @vars=(@vars,$name);
} #for

for($next_stack=$start_stack+1; $next_stack<=$end_stack+1; $next_stack++)
{
  while(!-e "outputs/$next_stack\_0000_elev.61") {sleep 300;} #sleep 5 min.
#  sleep 180; #wait a little longer to make sure outputs are written
  $current_stack=$next_stack-1;
  foreach $var (@vars)
  {
    if(!-e "outputs/$current_stack\_0000_$var") {next;}
    for($i=0;$i<$nproc;$i++)
    {
      $id=sprintf("%04d",$i);
      while(!-e "outputs/$next_stack\_$id\_$var") {sleep 180;} #sleep 3 min.
    } #for

    open(IN,">outputs/combine_output.in");
    print IN "$var\n $current_stack $current_stack\n $netcdf \n"; #3rd line inc=0: binary
    close(IN);
    system "cd outputs; $code >& combine.out";
    print "done combining $var in stack $current_stack...\n";
  } #foreach
} #for

#if(-e "outputs/$end_stack\_hvel.64") {system "cd outputs/; echo *_[0-9]???_*.6? | xargs rm -f";}
