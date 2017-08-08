set term pdf noenh size 1024,768
set outp "dmsc-cluster-sample-file.png"

set key t l
set title "DMSC cluster, Lustre stripe count 4, load 128 MByte chunk size 4 MByte"
set xla "MPI ranks"
set yla "MB/s" offset 1
p \
  "dmsc-cluster-sample-file-gzip" u 1:4 w p ti "GZIP" ls 100, \
  "dmsc-cluster-sample-file-no-compression" u 1:4 w p ti "no compression" ls 210
