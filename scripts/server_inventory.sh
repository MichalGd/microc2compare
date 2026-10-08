#!/usr/bin/env bash
set -uo pipefail
# Read-only inventory; unavailable commands are explicit and do not hide other facts.
probe() {
  printf '\n[%s]\n' "$*"
  if command -v "$1" >/dev/null 2>&1; then "$@" || printf 'FAILED: %s\n' "$*";
  else printf 'UNAVAILABLE: %s\n' "$1"; fi
}
probe cat /etc/os-release
probe uname -a
probe lscpu
probe lscpu -e=CPU,CORE,SOCKET,NODE,ONLINE
probe nproc
probe numactl --hardware
probe taskset -pc $$
probe free -b
probe cat /proc/self/cgroup
probe findmnt -t cgroup,cgroup2
probe lsblk -o NAME,SIZE,TYPE,ROTA,TRAN,FSTYPE,MOUNTPOINTS
probe df -hT
probe df -i
probe id
probe quota -s
probe uptime
probe ps -eo pid,user,comm,%cpu,%mem,rss --sort=-rss
probe vmstat 1 5
probe iostat -xz 1 5
# Environment enumeration can clean ~/.conda/environments.txt. Avoid it even
# when registration is disabled; inventory must leave existing home files alone.
probe conda --version
probe conda info --base
if [[ -r /proc/self/cgroup ]]; then
  while IFS=: read -r hierarchy controllers group; do
    if [[ "$hierarchy" == 0 ]]; then
      for file in cpu.max cpuset.cpus.effective memory.max memory.current; do
        [[ -r "/sys/fs/cgroup${group}/$file" ]] && probe cat "/sys/fs/cgroup${group}/$file"
      done
    fi
  done < /proc/self/cgroup
fi
printf '\nScheduler executables:\n'
for tool in sbatch squeue sinfo qsub qstat bsub bjobs; do command -v "$tool" || printf 'UNAVAILABLE: %s\n' "$tool"; done
