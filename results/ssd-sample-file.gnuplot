set term png noenh size 1024,768
set outp "ssd-sample-file.png"

set key t l
set title "DMSC workstation SSD, load 128 MByte chunk size 4 MByte"
set xla "MPI ranks"
set yla "MB/s" offset 1
p \
  "ssd-sample-file-gzip" u 1:4 w p ti "GZIP" ls 100, \
  "ssd-sample-file-no-compression" u 1:4 w p ti "no compression" ls 210
