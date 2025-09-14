package shared;

import java.io.PrintStream;

public class Timer {
	
	public Timer(){this(System.err, true);}
	
	public Timer(String s){
		this(System.err, true);
		if(outstream!=null){outstream.println(s);}
	}
	
	public Timer(PrintStream outstream_){
		this(outstream_, true);
	}
	
	public Timer(PrintStream outstream_, boolean addTab_){
		outstream=outstream_;
		addTab=addTab_;
		start();
	}
	
	public long start(String s){
		if(outstream!=null){outstream.println(s);}
		return start();
	}
	
	public long stopAndPrint(){
		long x=stop();
		if(outstream!=null){outstream.println(this);}
		return x;
	}
	
	public long stop(String s){
		long x=stop();
		if(addTab && s!=null && !s.endsWith("\t")){s=s+"\t";}
		if(outstream!=null){outstream.println(s+this);}
		return x;
	}
	
	public long stopAndStart(String s){
		long x=stop();
		if(addTab && s!=null && !s.endsWith("\t")){s=s+"\t";}
		if(outstream!=null){outstream.println(s+this);}
		start();
		return x;
	}
	
	public long start(){
		time1=time2=System.nanoTime();
		elapsed=0;
		return time1;
	}
	
	public long stop(){
		time2=System.nanoTime();
		elapsed=time2-time1;
		return time2;
	}
	
	@Override
	public String toString(){
//		new Exception().printStackTrace();
		return timeInSeconds(3)+" seconds.";
	}
	
	public String timeInSeconds(int decimals) {
		return Tools.format("%."+decimals+"f", timeInSeconds());
	}
	
	public double timeInSeconds() {
		return elapsed/1000000000d;
	}

	public long time1;
	public long time2;
	/** in nanos */
	public long elapsed;
	
	public PrintStream outstream=System.err;
	public boolean addTab=true;
}
