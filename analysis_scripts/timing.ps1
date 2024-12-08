$T1 = Measure-Command {./chiralplv2.exe EXAMPLE1thr.inp}
echo "1 thread"
echo ($T1.Milliseconds/1000)
$T4 = Measure-Command {./chiralplv2.exe EXAMPLE4thr.inp}
echo "4 threads"
echo ($T4.Milliseconds/1000)
$T8 = Measure-Command {./chiralplv2.exe EXAMPLE8thr.inp}
echo "8 threads"
echo ($T8.Milliseconds/1000)
$T16 = Measure-Command {./chiralplv2.exe EXAMPLE16thr.inp}
echo "16 threads"
echo ($T16.Milliseconds/1000)
$T32 = Measure-Command {./chiralplv2.exe EXAMPLE32thr.inp}
echo "32 threads"
echo ($T32.Milliseconds/1000)
$T64 = Measure-Command {./chiralplv2.exe EXAMPLE64thr.inp}
echo "64 threads"
echo ($T64.Milliseconds/1000)