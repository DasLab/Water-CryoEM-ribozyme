# Generation of Water network figures



```open ../../models/2.3A_SWIM.pdb
# style and coloring of the atomic models
set bgColor white
lighting soft
lighting shadows false
graphics silhouettes true
hide pseudobonds 
hide #9-10 cartoons
nucleotides #9,10 fill
style nucleic & #9-10 stick
color #9,10 #fff5ee
color #9,10:HOH #df0303
color #9,10:MG #00af2b
size #9,10:MG atomRadius 1
size #9,10:HOH atomRadius 0.8
style :HOH,MG sphere
show #9,10:HOH,MG a
# for panel A only the backbone was shown:
cartoon style ~:HOH,MG  xsect oval width 1.5 thick 1

# color density, zone near MG,
# color the cryo-EM densities
color zone #17 near #10:MG,HOH distance 2.5
vol splitbyzone #17
color #11.1,3 #ff000096
color #11.2 #00aa0096
color zone #16 near #9:MG,HOH distance 2.5
vol splitbyzone #16
color #12.1,3 #ff000096
color #12.2 #00aa0096
# color MD densities
vol #6 level 7.5 level 25 color #00aa0096 color #00aa00
vol #7 level 7.5 level 15 color #aa55ff96 color #aa55ff
vol #8 level 36.66666 level 55 color #ff000096 color #ff0000


# once found angle take images and save (B onwards)
hide all m
show #9,12 m,a
save Fig5B_2.2.png pixelSize 0.05 transparentBackground true
hide all m
show #10,11 m,a
save Fig5B_2.3.png pixelSize 0.05 transparentBackground true
hide all m
show #9,6,7,8 m,a
hide #9:HOH,MG a
save Fig5B_md.png pixelSize 0.05 transparentBackground true
save Fig5B.cxs
# A only see caroon of RNA
hide all m,a
show #9,12 m,c
cartoon style ~:HOH,MG  xsect oval width 1.5 thick 1
save Fig5A_2.2.png pixelSize 0.05 transparentBackground true
hide all m,a
show #10,11 m,c
save Fig5A_2.3.png pixelSize 0.05 transparentBackground true
hide all m,a
show #9,6,7,8 m,c
save Fig5A_md.png pixelSize 0.05 transparentBackground true
save Fig5A.cxs
```



## Generate Supplemental Movie 2

```
open Fig5B.cxs
hide all m
show #9,12 m,a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5B_2.2.mp4
hide all m
show #10,11 m,a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5B_2.3.mp4
hide all m
show #9,6,7,8 m,a
hide #9:HOH,MG a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5B_md.mp4

open Fig5C.cxs
hide all m
show #9,12 m,a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5C_2.2.mp4
hide all m
show #10,11 m,a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5C_2.3.mp4
hide all m
show #9,6,7,8 m,a
hide #9:HOH,MG a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5C_md.mp4

open Fig5D.cxs
hide all m
show #9,12 m,a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5D_2.2.mp4
hide all m
show #10,11 m,a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5D_2.3.mp4
hide all m
show #9,6,7,8 m,a
hide #9:HOH,MG a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5D_md.mp4

open Fig5E.cxs
hide all m
show #9,12 m,a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5E_2.2.mp4
hide all m
show #10,11 m,a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5E_2.3.mp4
hide all m
show #9,6,7,8 m,a
hide #9:HOH,MG a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5E_md.mp4

open Fig5F.cxs
hide all m
show #9,12 m,a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5F_2.2.mp4
hide all m
show #10,11 m,a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5F_2.3.mp4
hide all m
show #9,6,7,8 m,a
hide #9:HOH,MG a
movie record; turn x 120 200 wobble 200 wobbleaspect 0.3; wait 250; movie encode Fig5F_md.mp4
```

