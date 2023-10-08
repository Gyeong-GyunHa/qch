PROGRAM nondir_quadclustering

implicit none

integer:: i, n1, n2, size_list, ierr, io_num

integer, allocatable, dimension(:,:) :: inci_list
logical, allocatable, dimension(:,:) :: inci_mat
real, allocatable, dimension(:) :: quad_clu_coeff


! User should modify when changing input and output files
character (len=21) :: file_name
open(15,file='quadclustering_dist.txt')
file_name='sample_hypergraph.txt'


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

allocate (inci_list(size_list,2))
open(io_num,file=file_name)
do i=1,size_list
   read(io_num,*,iostat=ierr) inci_list(i,1),inci_list(i,2)
end do
close(io_num)

write(*,*) 'load done'
n1=maxval(inci_list(:,1))
n2=maxval(inci_list(:,2))

! transform the incidence list into incidence matrix
allocate (inci_mat(n1,n2))
inci_mat=.FALSE.
call inci_list2incmat( n1, n2, size_list, inci_list, inci_mat)

deallocate(inci_list)
allocate (quad_clu_coeff(n1))
write(*,*) 'start to compute the quad clustering coefficient'

! compute the quad clustering coefficient
call quad_clustering(n1, n2, quad_clu_coeff, inci_mat)

deallocate(inci_mat)

! export the distribution of the quad clustering coefficient
do i=1,n1
   write(15,*) quad_clu_coeff(i)
end do

deallocate(quad_clu_coeff)


close(15)
write(*,*) 'done'


END PROGRAM nondir_quadclustering



subroutine inci_list2incmat( n1, n2, size_list, inci_list, inci_mat)

   implicit none
   integer:: n1, n2, i, i1, i2, numofnodes1, size_list
   integer, allocatable, dimension(:) :: nodelist1
   logical, dimension(n1,n2) :: inci_mat
   integer, dimension(size_list,2) :: inci_list


   inci_mat=.FALSE.

   do i=1,n2
      numofnodes1=count(inci_list(:,2)==i)
      allocate(nodelist1(numofnodes1))
      i2=1
      do i1=findloc(inci_list(:,2),i,DIM=1),size_list
         if (inci_list(i1,2)==i) then
            nodelist1(i2)=inci_list(i1,1)
            i2=i2+1
         end if
         if (i2>numofnodes1) then
            exit
         end if
      end do
      !$OMP PARALLEL shared(inci_mat, nodelist1, i, numofnodes1) default(private)
         !$OMP do
            do i1=1,numofnodes1
               inci_mat(nodelist1(i1),i)=.TRUE.
            end do
         !$OMP end do
      !$OMP END PARALLEL
      deallocate(nodelist1)
   end do

end subroutine inci_list2incmat


subroutine quad_clustering(n1, n2, clu_coeff, inci_mat)

   implicit none
   integer:: n1, n2, i, denominator, numerator, numofedges
   integer:: a, b, j, i1, i2, i3, i4
   integer, allocatable, dimension(:) :: edgelist1
   real, dimension(n1) :: clu_coeff
   logical, dimension(n1,n2) :: inci_mat


   clu_coeff=0.0

   !$OMP PARALLEL shared(n1, n2, clu_coeff, inci_mat) default(private)
      !$OMP do
         do i=1,n1
            denominator=0
            numerator=0
            numofedges=count(inci_mat(i,:).eqv..TRUE.)
            if (numofedges<2) then
               clu_coeff(i)=0.0
               cycle
            end if
            allocate(edgelist1(numofedges))
            edgelist1=0
            i2=1
            do a=1,n2
               if (inci_mat(i,a).eqv..TRUE.) then
                  edgelist1(i2)=a
                  i2=i2+1
               end if
               if (i2>numofedges) then
                  exit
               end if
            end do
            do a=1,numofedges-1
               if (count(inci_mat(:,edgelist1(a)).eqv..TRUE.)==1) cycle
               i1=findloc(inci_mat(:,edgelist1(a)),.TRUE.,DIM=1)
               i2=findloc(inci_mat(:,edgelist1(a)),.TRUE.,DIM=1,BACK=.TRUE.)
               do b=a+1,numofedges
                  if (count(inci_mat(:,edgelist1(b)).eqv..TRUE.)==1) cycle
                  i3=findloc(inci_mat(:,edgelist1(b)),.TRUE.,DIM=1)
                  i4=findloc(inci_mat(:,edgelist1(b)),.TRUE.,DIM=1,BACK=.TRUE.)
                  denominator=denominator+minval((/count(inci_mat(:,edgelist1(a)).eqv..TRUE.)-1, &
count(inci_mat(:,edgelist1(b)).eqv..TRUE.)-1/))
                  do j=maxval((/i1, i3/)),minval((/i2, i4/))
                     if ((inci_mat(j,edgelist1(a)).eqv..TRUE.) .and. (inci_mat(j,edgelist1(b)).eqv..TRUE.) .and. j/=i) then
                        numerator=numerator+1
                     end if
                  end do
               end do
            end do
            if (denominator==0) then
               clu_coeff(i)=0
            else
               clu_coeff(i)=(numerator*1.0)/(denominator*1.0)
            end if
            deallocate(edgelist1)
         end do
      !$OMP end do
   !$OMP END PARALLEL

end subroutine quad_clustering