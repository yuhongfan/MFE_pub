Module blkdata
       implicit none
       save
!
       integer :: iflg_blk, np_blk, mp_blk, n_blk, m_blk
       integer :: idimy_blk, nww
       integer :: kk,ll,ierr_blk
       real(8), allocatable :: an(:),bn(:),cn(:)
       real(8), allocatable :: am(:),bm(:),cm(:)
       real(8), allocatable :: yy(:,:),ww(:)
!
       integer :: nrlpts
       integer :: kint
       real(8) :: phk
       real(8), allocatable :: rldat(:)
       real(8), allocatable ::  wsave(:)
       real(8), allocatable ::  qq(:, :, :)
!
End Module blkdata
