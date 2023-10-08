PROGRAM dir_quadclustering

implicit none

integer:: i, n1, n2, size_list, ierr, io_num
integer, allocatable, dimension(:,:) :: inci_list
logical, allocatable, dimension(:,:) :: inci_mat_forward, inci_mat_backward
real, allocatable, dimension(:) :: quad_clu_coeff


! User should modify when changing input and output files
character (len=21) :: file_name
open(15,file='dirquadclustering_dist.txt')
file_name='dir_hyperg_sample.txt'

! examine the size of input file
io_num=20
open(io_num,file=file_name)
i=1
do 
   read(io_num,*,iostat=ierr)
   if (ierr<0) exit
   i=i+1
end do
close(io_num)
size_list=i-1

allocate (inci_list(size_list,3))
open(io_num,file=file_name)
do i=1,size_list
   read(io_num,*,iostat=ierr) inci_list(i,1), inci_list(i,2), inci_list(i,3)
end do
close(io_num)

write(*,*) 'load done'

n1=maxval(inci_list(:,1))
n2=maxval(inci_list(:,2))

allocate (inci_mat_forward(n1,n2))
allocate (inci_mat_backward(n1,n2))
call inci_list2dir_incimat( n1, n2, size_list, inci_list, inci_mat_forward, inci_mat_backward)

deallocate(inci_list)
allocate (quad_clu_coeff(n1))
write(*,*) 'start to compute the directed quad clustering coefficient'

! compute the quad clustering coefficient
call dir_quad_clustering(n1, n2, quad_clu_coeff, inci_mat_forward, inci_mat_backward)

deallocate(inci_mat_forward)
deallocate(inci_mat_backward)

! export the distribution of the quad clustering coefficient
do i=1,n1
   write(15,*) quad_clu_coeff(i)
end do

deallocate(quad_clu_coeff)

close(15)
write(*,*) 'done'


END PROGRAM dir_quadclustering




subroutine inci_list2dir_incimat( n1, n2, size_list, inci_list, inci_mat_forward, inci_mat_backward)

   implicit none
   integer:: n1, n2, i, size_list
   logical, dimension(n1,n2) :: inci_mat_forward, inci_mat_backward
   integer, dimension(size_list,3) :: inci_list


   inci_mat_forward=.FALSE.
   inci_mat_backward=.FALSE.
   !$OMP PARALLEL shared(size_list, inci_mat_forward, inci_mat_backward, inci_list)
      !$OMP do
         do i=1,size_list
            if (inci_list(i,3)==1) then
               inci_mat_forward(inci_list(i,1),inci_list(i,2))=.TRUE.
            else
               inci_mat_backward(inci_list(i,1),inci_list(i,2))=.TRUE.
            end if
         end do
      !$OMP end do
   !$OMP END PARALLEL

end subroutine inci_list2dir_incimat








subroutine dir_quad_clustering(n1, n2, clu_coeff, inci_mat_forward, inci_mat_backward)

   implicit none
   integer:: n1, n2, numofnodes1, numofnodes2, numofnodes3, numofnodes4
   integer:: i, i1, i2, j, denominator, numerator
   logical :: bool1, bool2
   logical, dimension(n1,n2) :: inci_mat_forward, inci_mat_backward
   real, dimension(n1) :: clu_coeff

   clu_coeff=0.0


   numofnodes1=0 ! alpha j forward
   numofnodes2=0 ! alpha j backward
   numofnodes3=0 ! beta j forward
   numofnodes4=0 ! beta j baclward

   do i=1,n1 ! node i
      denominator=0
      numerator=0
      do i1=1,n2-1 ! edge alpha
         if ((inci_mat_forward(i,i1).eqv..FALSE.) .and. (inci_mat_backward(i,i1).eqv..FALSE.)) cycle
         do i2=i1+1,n2 ! edge beta
            if ((inci_mat_forward(i,i2).eqv..FALSE.) .and. (inci_mat_backward(i,i2).eqv..FALSE.)) cycle    
            do j=1,n1 ! node j
               bool1=.FALSE.
               bool2=.FALSE.
               if ((inci_mat_forward(j,i1).eqv..FALSE.) .and. (inci_mat_backward(j,i1).eqv..FALSE.)) bool1=.TRUE.
               if ((inci_mat_forward(j,i2).eqv..FALSE.) .and. (inci_mat_backward(j,i2).eqv..FALSE.)) bool2=.TRUE.
               if (i==j) cycle
               if ((bool1.eqv..TRUE.) .and. (bool2.eqv..TRUE.)) then
                  cycle
               else if ((bool1.eqv..FALSE.) .and. (bool2.eqv..TRUE.)) then ! j,i1 alpha 에 대하여 
                  numofnodes1=numofnodes1+(1*count((/inci_mat_forward(j,i1)/).eqv..TRUE.))
                  numofnodes2=numofnodes2+(1*count((/inci_mat_backward(j,i1)/).eqv..TRUE.))
               else if ((bool1.eqv..TRUE.) .and. (bool2.eqv..FALSE.)) then ! j,i2 beta 에 대하여 
                  numofnodes3=numofnodes3+(1*count((/inci_mat_forward(j,i2)/).eqv..TRUE.))
                  numofnodes4=numofnodes4+(1*count((/inci_mat_backward(j,i2)/).eqv..TRUE.))
               else if ((bool1.eqv..FALSE.) .and. (bool2.eqv..FALSE.)) then ! j,i1 alpha, i2 beta 모두에 대하여
                  numofnodes1=numofnodes1+(1*count((/inci_mat_forward(j,i1)/).eqv..TRUE.))
                  numofnodes2=numofnodes2+(1*count((/inci_mat_backward(j,i1)/).eqv..TRUE.))
                  numofnodes3=numofnodes3+(1*count((/inci_mat_forward(j,i2)/).eqv..TRUE.))
                  numofnodes4=numofnodes4+(1*count((/inci_mat_backward(j,i2)/).eqv..TRUE.))
                  numerator=numerator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)&
*count((/inci_mat_forward(j,i1), inci_mat_backward(j,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(j,i2), inci_mat_backward(j,i2)/).eqv..TRUE.))
               end if
            end do
            if ((numofnodes1==numofnodes2) .and. (numofnodes3==numofnodes4)) then
               denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*(4*minval((/numofnodes1,numofnodes3/))))
            else if ((numofnodes1==numofnodes2) .and. (numofnodes3/=numofnodes4)) then
               if (minval((/numofnodes1,numofnodes3,numofnodes4/))==numofnodes1) then
                  denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*(4*numofnodes1))
               else if ((numofnodes1/=numofnodes3) .and. (numofnodes1/=numofnodes4)) then
                  if (minval((/numofnodes1,numofnodes3,numofnodes4/))==numofnodes3) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((2*numofnodes3)&
+(2*minval((/numofnodes1,numofnodes4/)))))
                  else if (minval((/numofnodes1,numofnodes3,numofnodes4/))==numofnodes4) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((2*numofnodes4)&
+(2*minval((/numofnodes1,numofnodes3/)))))
                  end if
               else
                  denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((2*minval((/numofnodes1,numofnodes3,numofnodes4/)))&
+(2*maxval((/numofnodes1,numofnodes3,numofnodes4/)))))
               end if
            else if ((numofnodes1/=numofnodes2) .and. (numofnodes3==numofnodes4)) then
               if (minval((/numofnodes3,numofnodes1,numofnodes2/))==numofnodes3) then
                  denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*(4*numofnodes3))
               else if ((numofnodes3/=numofnodes1) .and. (numofnodes3/=numofnodes2)) then
                  if (minval((/numofnodes3,numofnodes1,numofnodes2/))==numofnodes1) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((2*numofnodes1)&
+(2*minval((/numofnodes3,numofnodes2/)))))
                  else if (minval((/numofnodes3,numofnodes1,numofnodes2/))==numofnodes2) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((2*numofnodes2)&
+(2*minval((/numofnodes3,numofnodes1/)))))
                  end if
               else
                  denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((2*minval((/numofnodes3,numofnodes1,numofnodes2/)))&
+(2*maxval((/numofnodes3,numofnodes1,numofnodes2/)))))
               end if
            else if (((numofnodes1/=numofnodes2) .and. (numofnodes3/=numofnodes4)) .and. ((numofnodes1==numofnodes3) .or. &
(numofnodes1==numofnodes4) .or. (numofnodes2==numofnodes3) .or. (numofnodes2==numofnodes4))) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               if ((minval((/numofnodes1,numofnodes2/))==minval((/numofnodes3,numofnodes4/))) .and. &
(maxval((/numofnodes1,numofnodes2/))==maxval((/numofnodes3,numofnodes4/)))) then
                  denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*minval((/numofnodes1,numofnodes2/)))&
+(1*(maxval((/numofnodes1,numofnodes2/))-minval((/numofnodes1,numofnodes2/))))))
               else if ((minval((/numofnodes1,numofnodes2/))==minval((/numofnodes3,numofnodes4/))) .and. &
(maxval((/numofnodes1,numofnodes2/))/=maxval((/numofnodes3,numofnodes4/)))) then
                  if (numofnodes1==numofnodes3) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*numofnodes1)&
+(1*(minval((/numofnodes2,numofnodes4/))-numofnodes1))))
                  else if (numofnodes2==numofnodes3) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*numofnodes2)&
+(1*(minval((/numofnodes1,numofnodes4/))-numofnodes2))))
                  else if (numofnodes1==numofnodes4) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*numofnodes1)&
+(1*(minval((/numofnodes2,numofnodes3/))-numofnodes1))))
                  else if (numofnodes2==numofnodes4) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*numofnodes2)&
+(1*(minval((/numofnodes1,numofnodes3/))-numofnodes2))))
                  end if
               else if (maxval((/numofnodes1,numofnodes2/))==minval((/numofnodes3,numofnodes4/))) then
                  if (numofnodes1==numofnodes3) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*numofnodes2)&
+(2*(numofnodes1-numofnodes2))))
                  else if (numofnodes2==numofnodes3) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*numofnodes1)&
+(2*(numofnodes2-numofnodes1))))
                  else if (numofnodes1==numofnodes4) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*numofnodes2)&
+(2*(numofnodes1-numofnodes2))))
                  else if (numofnodes2==numofnodes4) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*numofnodes1)&
+(2*(numofnodes2-numofnodes1))))
                  end if
               else if (minval((/numofnodes1,numofnodes2/))==maxval((/numofnodes3,numofnodes4/))) then
                  if (numofnodes1==numofnodes3) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*numofnodes4)&
+(2*(numofnodes3-numofnodes4))))
                  else if (numofnodes2==numofnodes3) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*numofnodes4)&
+(2*(numofnodes3-numofnodes4))))
                  else if (numofnodes1==numofnodes4) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*numofnodes3)&
+(2*(numofnodes4-numofnodes3))))
                  else if (numofnodes2==numofnodes4) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*numofnodes3)&
+(2*(numofnodes4-numofnodes3))))
                  end if
               else if ((minval((/numofnodes1,numofnodes2/))/=minval((/numofnodes3,numofnodes4/))) .and. &
(maxval((/numofnodes1,numofnodes2/))==maxval((/numofnodes3,numofnodes4/)))) then
                  if (numofnodes1==numofnodes3) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*minval((/numofnodes2,numofnodes4/)))&
+(2*(maxval((/numofnodes2,numofnodes4/))-minval((/numofnodes2,numofnodes4/))))&
+(1*(numofnodes1-maxval((/numofnodes2,numofnodes4/))))))
                  else if (numofnodes2==numofnodes3) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*minval((/numofnodes1,numofnodes4/)))&
+(2*(maxval((/numofnodes1,numofnodes4/))-minval((/numofnodes1,numofnodes4/))))&
+(1*(numofnodes2-maxval((/numofnodes1,numofnodes4/))))))
                  else if (numofnodes1==numofnodes4) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*minval((/numofnodes2,numofnodes3/)))&
+(2*(maxval((/numofnodes2,numofnodes3/))-minval((/numofnodes2,numofnodes3/))))&
+(1*(numofnodes1-maxval((/numofnodes2,numofnodes3/))))))
                  else if (numofnodes2==numofnodes4) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((4*minval((/numofnodes1,numofnodes3/)))&
+(2*(maxval((/numofnodes1,numofnodes3/))-minval((/numofnodes1,numofnodes3/))))&
+(1*(numofnodes2-maxval((/numofnodes1,numofnodes3/))))))
                  end if
               end if
            else if ((numofnodes1/=numofnodes2) .and. (numofnodes3/=numofnodes4) .and. (numofnodes1/=numofnodes3) .and. &
(numofnodes1/=numofnodes4) .and. (numofnodes2/=numofnodes3) .and. (numofnodes2/=numofnodes4)) then
               if (maxval((/numofnodes1,numofnodes2/))<minval((/numofnodes3,numofnodes4/))) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((2*minval((/numofnodes1,numofnodes2/)))&
+(2*maxval((/numofnodes1,numofnodes2/)))))
               else if (minval((/numofnodes1,numofnodes2/))<minval((/numofnodes3,numofnodes4/))) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((2*minval((/numofnodes1,numofnodes2/)))&
+(1*minval((/maxval((/numofnodes1,numofnodes2/)),minval((/numofnodes3,numofnodes4/))/)))&
+(1*maxval((/maxval((/numofnodes1,numofnodes2/)),minval((/numofnodes3,numofnodes4/))/)))))
               else if (maxval((/numofnodes3,numofnodes4/))<minval((/numofnodes1,numofnodes2/))) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((2*minval((/numofnodes3,numofnodes4/)))&
+(2*maxval((/numofnodes3,numofnodes4/)))))
               else if (minval((/numofnodes3,numofnodes4/))<minval((/numofnodes1,numofnodes2/))) then
                     denominator=denominator+(1*count((/inci_mat_forward(i,i1), inci_mat_backward(i,i1)/).eqv..TRUE.)&
*count((/inci_mat_forward(i,i2), inci_mat_backward(i,i2)/).eqv..TRUE.)*((2*minval((/numofnodes3,numofnodes4/)))&
+(1*minval((/maxval((/numofnodes3,numofnodes4/)),minval((/numofnodes1,numofnodes2/))/)))&
+(1*maxval((/maxval((/numofnodes3,numofnodes4/)),minval((/numofnodes1,numofnodes2/))/)))))
               end if
            end if
            numofnodes1=0
            numofnodes2=0
            numofnodes3=0
            numofnodes4=0
         end do
      end do
      if (denominator==0) then
         clu_coeff(i)=0.0
      else
         clu_coeff(i)=(numerator*1.0)/(denominator*1.0)
      end if
   end do



end subroutine dir_quad_clustering