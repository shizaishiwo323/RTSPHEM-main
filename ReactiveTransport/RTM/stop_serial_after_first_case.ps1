param(
    [Parameter(Mandatory)] [string] $StdoutPath,
    [Parameter(Mandatory)] [int] $ParentPid,
    [Parameter(Mandatory)] [int] $ChildPid,
    [Parameter(Mandatory)] [string] $AuditLog
)

while (Get-Process -Id $ChildPid -ErrorAction SilentlyContinue) {
    if ((Get-Content -Raw $StdoutPath -ErrorAction SilentlyContinue) -match '\[2/4\]\s+M14_25') {
        "$(Get-Date -Format o) detected serial case 2; stopping PIDs $ChildPid,$ParentPid" |
            Set-Content -LiteralPath $AuditLog
        Stop-Process -Id $ChildPid -Force -ErrorAction SilentlyContinue
        Stop-Process -Id $ParentPid -Force -ErrorAction SilentlyContinue
        exit 0
    }
    Start-Sleep -Milliseconds 200
}

"$(Get-Date -Format o) main child exited before serial case 2" |
    Set-Content -LiteralPath $AuditLog
