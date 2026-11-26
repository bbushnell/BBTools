package jgi;

import java.io.File;
import java.util.ArrayList;

import parse.Parse;
import shared.Tools;

/**
 * Runs AssemblyStats2 sequentially for an ordered list of assembly inputs.
 * Readable positional arguments without '=' are accepted wherever they occur.
 * Named in/ref values naming an existing path stay literal; other values are
 * split on commas. Named inputs without a value add no inputs. Remaining
 * arguments are forwarded to every invocation.
 * Supplies header=t, showspeed=f, addname=t and k=0 before forwarded options;
 * later caller options can override these defaults. For each input after the
 * first, changes the wrapper's header default to header=f, calls System.gc(),
 * waits with a 100-millisecond timeout, and yields. Collection and a fixed pause
 * are not guaranteed, and there is no input-size threshold. Interruption of
 * the wait is ignored.
 * Constructs a new AssemblyStats2 for each input. After constructing each
 * instance after the first, sets its shared overwrite=false and append=true
 * before processing; other AssemblyStats2 static state is not reset here.
 *
 * @author Brian Bushnell
 * @date Apr 17, 2013
 */
public class AssemblyStatsWrapper {
	
	/**
	 * Collects inputs and invokes AssemblyStats2 once for each, in encounter order.
	 * @param args Input paths, named in/ref lists, and options forwarded to each input
	 */
	public static void main(String[] args){
		ArrayList<String> alist=new ArrayList<String>();
		ArrayList<String> ilist=new ArrayList<String>();
		
		alist.add("");
		alist.add("header=t");
		alist.add("showspeed=f");
		alist.add("addname=t");
		alist.add("k=0");
		
		for(String arg : args){
			if(!arg.contains("=") && Tools.canRead(arg)){
				ilist.add("in="+arg);
			}else{
				//JGI-058: later '=' characters belong to the filename, not the option key.
				String[] split=Parse.splitOnFirst(arg, '=');
				if(split[0].equalsIgnoreCase("in") || split[0].equalsIgnoreCase("ref")){
					if(split.length>1){
						if(new File(split[1]).exists()){
							ilist.add(arg);
						}else{
							String[] split2=split[1].split(",");
							for(String s : split2){
								ilist.add("in="+s);
							}
						}
					}
				}else{
					alist.add(arg);
				}
			}
		}
		
		String[] args2=alist.toArray(new String[alist.size()]);
		for(int i=0; i<ilist.size(); i++){
			String s=ilist.get(i);
//			System.err.println("Processing "+s);
			args2[0]=s;
			if(i>0){
				args2[1]="header=f";
//				AssemblyStats2.reset();
				System.gc();
				synchronized(AssemblyStatsWrapper.class){
					try {
						AssemblyStatsWrapper.class.wait(100);
					} catch (InterruptedException e) {}
				}
				Thread.yield();
			}
			AssemblyStats2 as2=new AssemblyStats2(args2);
			if(i>0){
				AssemblyStats2.overwrite=false;
				AssemblyStats2.append=true;
			}
			as2.process();
		}
		
	}
	
}
