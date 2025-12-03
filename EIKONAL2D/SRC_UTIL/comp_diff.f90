program diff

real*4 :: tab1(501,501),tab2(501,501),tab(501,501)

open(7,file='TT.bin',access='direct',status='old',recl=501*501*4)
open(8,file='TT_ref.bin',access='direct',status='old',recl=501*501*4)

read(7,rec=1) tab1
read(8,rec=1) tab2

close(7)
close(8)

tab(:,:)=tab1(:,:)-tab2(:,:)

open(9,file='TT_diff.bin',access='direct',status='unknown',recl=501*501*4)

write(9,rec=1) tab

close(9)

end program diff
