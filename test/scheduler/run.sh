#!/bin/sh
#Run the scheduler tests (see docs/94_scheduler.md). It starts the scheduler in
#BINDIR with a private home folder and a free port, and stops it at the end.
#Usage: run.sh SRCDIR BINDIR PROBE [lws]
src=$1
bin=$2
probe=$3
lws=$4
dir=$(cd "$(dirname "$0")" && pwd)
if [ ! -x "$bin/scheduler" ];then
	echo "$bin/scheduler is not built. Skip."
	exit 0
fi
home=$(mktemp -d "${TMPDIR:-/tmp}/maos_scheduler_test.XXXXXX") || exit 1
#a port that is free together with port+100 used by libwebsockets
port=$(python3 -c '
import random, socket
while True:
	p=random.randint(20000, 30000)
	try:
		for q in (p, p+100):
			s=socket.socket(); s.bind(("127.0.0.1", q)); s.close()
		print(p); break
	except OSError:
		pass')
mkdir -p "$home/.aos"
echo "$port" > "$home/.aos/port"
(cd "$home" && HOME="$home" exec "$bin/scheduler" > "$home/scheduler.log" 2>&1) &
pid_sch=$!
for i in 1 2 3 4 5 6 7 8 9 10;do
	ls "$home"/.aos/tmp/*/scheduler >/dev/null 2>&1 && break
	sleep 1
done
ans=0
for test in "test_commands.py $home" "test_auth.py $home $bin $probe $src" "test_websocket.py $home";do
	echo "== $test"
	python3 -B "$dir"/$test || ans=1
done
if [ "$lws" = 1 ];then
	echo "== test_websocket.py $home lws"
	python3 -B "$dir/test_websocket.py" "$home" lws || ans=1
fi
kill "$pid_sch"
wait "$pid_sch" 2>/dev/null
if [ $ans != 0 ];then
	echo "Scheduler log: $home/scheduler.log"
else
	rm -rf "${home:?}"
fi
exit $ans
