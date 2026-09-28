! State shared with the B1 investigation hooks (rt_hooks.f90). Not for upstream.
module rt_state
    use common_types, only: t_model_arguments
    implicit none
    type(t_model_arguments), save :: sv_model
    double precision, save :: sv_mudisk = 0.d0, sv_risco = 0.d0, sv_mueff = 0.d0
    integer, save :: sv_nf = 0, sv_xe = 0, sv_ne = 0
    double precision, allocatable, save :: sv_fi(:)
    ! copy of the kernel produced by rtrans (nlp = 1, me = 1): (ne, nf, xe)
    complex, allocatable, save :: sv_w0(:,:,:), sv_w1(:,:,:), sv_w2(:,:,:), sv_w3(:,:,:)
    ! externally supplied kernel that replaces rtrans's own
    logical, save :: ov_active = .false.
    complex, allocatable, save :: ov_w0(:,:,:), ov_w1(:,:,:)
end module rt_state
