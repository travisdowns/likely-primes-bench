; CONST SEGMENT READONLY ALIGN(32) 'CONST'
%use altreg

SECTION .rodata align=32
; We are only using 30 primes, so the last 2 are left as zero
allones:
db 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
db 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0

; The first 30 prime numbers (minus 1)
asmprimes:
db   2,   4,   6,  10,  12,  16,  18,  22
db  28,  30,  36,  40,  42,  46,  52,  58
db  60,  66,  70,  72,  78,  82,  88,  96
db 100, 102, 106, 108, 112, 126,   0,   0

; CONST ENDS

; .code

%macro pushall 0
push r0
push r1
push r2
push r3
push r4
push r5
push r6
push r7
push r8
push r9
push r10
push r12
push r13
push r14
push r15
%endmacro

%macro popall 0
pop r15
pop r14
pop r13
pop r12
pop r10
pop r9
pop r8
pop r7
pop r6
pop r5
pop r4
pop r3
pop r2
pop r1
pop r0
%endmacro




section .text

; Call with:
; extern "C" ULONGLONG ProcessA(ULONGLONG start, ULONGLONG end, const BYTE *state);

; ProcessA proc
GLOBAL ProcessA_asm:function

ProcessA_asm:
mov r15, rcx  ; save the callback struct address in r15

; map from the Linux registers to the windows ones that are assumed below
mov r8, rdx
mov r9, rcx
mov rcx, rdi
mov rdx, rsi




; Windows calling convention:

; rcx - start
; rdx - end
; r8 - current state
; r9, r10, r11 - scratch
; rax - return value

vmovups ymm0, [allones]     ; Set all bytes to 1
vmovdqu ymm1, [asmprimes]   ; Load primes
vmovdqu ymm2, [r8]          ; Load the current state

xor eax, eax                            ; Returns # of prime-ish found

xor r9d, r9d                            ; How many steps to next prime-ish?
vpsubb ymm2, ymm2, ymm0                 ; Decrement all the bytes by 1
jmp .top

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; This loop is the bottleneck.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
align 16

.nottop:

add r9d, 1
vpminub ymm2, ymm2, ymm1                ; Adjust any underflows
vpsubb ymm2, ymm2, ymm0                 ; Decrement all the bytes by 1

.top:
vpmovmskb r10d, ymm2                    ; Check the upper bits for underflows

test r10d, r10d
jnz .nottop                              ; No underflows means prime-ish
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of bottleneck loop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; r9 contains the count of iterations (0-97 when using 30 primes)

lea r9, [r9 + r9 + 2]                   ; Account for even numbers and starting at 0

; Decrement the current value to the next 'likely' prime.
; This may require sbb for values greater than 64 bits.
sub rcx, r9

; At this point, rcx contains the likely prime.

cmp rcx, rdx                            ; if (start < end) goto done;
jb .done

;pushall
;mov rdi, [r15 + 8]
;call [r15]
;popall

; The vpsubb is done here in the (mostly vain) hope that it can be scheduled
; to run concurrently with other instructions.

xor r9d, r9d                            ; Reset for next loop
vpsubb ymm2, ymm2, ymm0                 ; Decrement all the bytes by 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; This is a placeholder for a more involved routine.
;
; There's no simple way to provide a runable version
; of it, since it makes use of a lookup table that
; requires libraries to create and is too big (8gig)
; to upload.
;
; OTOH, if the loop above can be optimized to the
; point where *this* code becomes the bottleneck, I'll
; see what I can do.
;
; FYI it makes no use of any xmm or ymm register, nor
; does it 'call' any other routines (it's just math).
; It does make use of most GP registers.

inc rax ; Count how many we find

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

jmp .top

.done:

; rax already contains the return value.
; All the registers used are 'scratch', so didn't get saved and
; don't need to be restored.
; All that leaves is clearing the ymm registers.
vzeroall

ret

GLOBAL ProcessA2_asm:function

ProcessA2_asm:
mov r15, rcx  ; save the callback in r15

; map from the Linux registers to the windows ones that are assumed below
mov r8, rdx
mov r9, rcx
mov rcx, rdi
mov rdx, rsi




; Windows calling convention:

; rcx - start
; rdx - end
; r8 - current state
; r9, r10, r11 - scratch
; rax - return value

vmovups ymm0, [allones]     ; Set all bytes to 1
vmovdqu ymm1, [asmprimes]   ; Load primes
vmovdqu ymm2, [r8]          ; Load the current state

xor eax, eax                            ; Returns # of prime-ish found

xor r9d, r9d                            ; How many steps to next prime-ish?
vpsubb ymm2, ymm2, ymm0                 ; Decrement all the bytes by 1

sub rcx, rdx
shr rcx, 1

jmp .top

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; This loop is the bottleneck.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
align 32

.top:

vpminub ymm2, ymm2, ymm1                ; Adjust any underflows
vpsubb ymm2, ymm2, ymm0                 ; Decrement all the bytes by 1

vpmovmskb r10d, ymm2                    ; Check the upper bits for underflows

test r10d, r10d

setz r14l
add rax, r14

dec rcx
jnz .top



;pushall
;call r15
;popall


; rax already contains the return value.
; All the registers used are 'scratch', so didn't get saved and
; don't need to be restored.
; All that leaves is clearing the ymm registers.
vzeroupper

ret


; ProcessA endp

;end

EXTERN BYTE_BITMAPS256
EXTERN BYTE_OFFSETS256
EXTERN OFFSET_PERIODS256

EXTERN BYTE_BITMAPS512
EXTERN BYTE_OFFSETS512
EXTERN OFFSET_PERIODS512

GLOBAL kernel64_asm:function

kernel64_asm:
push   rbp
mov    rbp,rsp
push   r15
push   r14
push   r13
push   r12
push   r10
push   rbx
lea    r13,[rsi+0x400]
mov    r12,rsi
xor    ebx,ebx
.top:
mov    r15d, BYTE_BITMAPS256
xor    ecx,ecx
xor    esi,esi
xor    r8d,r8d
xor    r9d,r9d
xor    r10d,r10d
xor    r11d,r11d
.inner:
movzx  r14d,BYTE [rdi+rcx*1]
mov    eax,ecx
shl    rax,0x7
movzx  edx,BYTE [r14 + rax + BYTE_OFFSETS256]

or     r8, [r15 + rdx]
or     r9, [r15 + rdx + 8]
or     r10,[r15 + rdx + 16]
or     r11,[r15 + rdx + 24]

lea    r14d,[r14+0x1]
cmp    r14l,BYTE [OFFSET_PERIODS256 + rcx]
cmove  r14d,ebx
add    r15,0x9e
mov    BYTE [rdi+rcx*1],r14l
add    rcx,0x1
cmp    rcx,0x1e
jne    .inner
add    r12,0x20
vmovq  xmm3,r8
vpinsrq xmm0,xmm3,r9,0x1
vmovq  xmm2,r10
vpinsrq xmm1,xmm2,r11,0x1
vinserti128 ymm0,ymm0,xmm1,0x1
vmovups [r12-0x20],ymm0

cmp    r13,r12
jne    .top

vzeroupper
pop    rbx
pop    r10
pop    r12
pop    r13
pop    r14
pop    r15
mov    rsp, rbp
pop    rbp
ret

GLOBAL kernel256_asm:function

kernel256_asm:
    push   rbp
    mov    rbp,rsp
    push   r15
    push   r14
    push   r13
    push   r12
    push   r10
    push   rbx
    lea    r13,[rsi+0x400]
    mov    r12,rsi
    xor    ebx,ebx
.top:
    mov    r15d, BYTE_BITMAPS256
    xor    ecx,ecx
    xor    esi,esi
    vpxor  xmm0, xmm0
.inner:
    movzx  r14d,BYTE [rdi + rcx]
    mov    eax,ecx
    shl    rax,0x7
    movzx  edx,BYTE [r14 + rax + BYTE_OFFSETS256]

    vpor   ymm0, ymm0, [r15 + rdx]

    lea    r14d,[r14+0x1]
    cmp    r14l,BYTE [OFFSET_PERIODS256 + rcx]
    cmove  r14d,ebx
    add    r15,0x9e
    mov    BYTE [rdi+rcx*1],r14l
    add    rcx,1
    cmp    rcx,30
    jne    .inner

    vmovups [r12], ymm0
    add    r12,0x20

    cmp    r13,r12
    jne    .top

    vzeroupper
    pop    rbx
    pop    r10
    pop    r12
    pop    r13
    pop    r14
    pop    r15
    mov    rsp, rbp
    pop    rbp
    ret

GLOBAL kernel512_asm:function

kernel512_asm:
    push   rbp
    mov    rbp,rsp
    push   r15
    push   r14
    push   r13
    push   r12
    push   r10
    push   rbx
    lea    r13,[rsi+0x800]
    mov    r12,rsi
    xor    ebx,ebx
.top:
    mov    r15d, BYTE_BITMAPS512
    xor    ecx,ecx
    xor    esi,esi
    vpxor  xmm0, xmm0
    vpxor  xmm1, xmm1
.inner:
    movzx  r14d,BYTE [rdi + rcx]
    mov    eax,ecx
    shl    rax,0x7
    movzx  edx,BYTE [r14 + rax + BYTE_OFFSETS512]

    vpor   ymm0, ymm0, [r15 + rdx]
    vpor   ymm1, ymm1, [r15 + rdx + 32]

    lea    r14d,[r14+0x1]
    cmp    r14l,BYTE [OFFSET_PERIODS512 + rcx]
    cmove  r14d,ebx
    add    r15,190
    mov    BYTE [rdi+rcx*1],r14l
    add    rcx,1
    cmp    rcx,30
    jne    .inner

    vmovups [r12],    ymm0
    vmovups [r12+32], ymm1
    add    r12,64

    cmp    r13,r12
    jne    .top

    vzeroupper
    pop    rbx
    pop    r10
    pop    r12
    pop    r13
    pop    r14
    pop    r15
    mov    rsp, rbp
    pop    rbp
    ret





