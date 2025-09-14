package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Parse;
import shared.PreParser;
import shared.Shared;
import shared.Timer;

/**
 * Renames files based on their headers
 * @author Brian Bushnell
 * @date May 19, 2016
 *
 */
public class RenameByHeader {
	
	public static void main(String[] args){
		Timer t=new Timer();
		RenameByHeader x=new RenameByHeader(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public RenameByHeader(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=false;
		ReadWrite.USE_UNPIGZ=false;
		ReadWrite.USE_UNBGZIP=false;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			File f=(b==null ? new File(arg) : null);

			if(f!=null && f.exists()){
				if(f.isDirectory()){
					for(File f2 : f.listFiles()){
						String name=f2.getAbsolutePath();
						if(f2.isFile() && FileFormat.hasFastqOrFastqExtension(name)){
							list.add(name);
						}
					}
				}else{
					list.add(f.getAbsolutePath());
				}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
	}
	
	void process(Timer t){
		for(String s : list){
			processFile(s);
		}
	}
	
	void processFile(String path){
		TextFile tf=new TextFile(path);
		String line=tf.nextLine();
		tf.close();
		if(line==null){return;}
		
		StringBuilder sb=new StringBuilder();
		File f=new File(path);
		String dir=f.getParent();
		if(dir!=null){sb.append(dir).append('/');}
		try {
			String[] split=line.substring(1).replace(",", "").split(" ");
			sb.append(split[1]);
			sb.append('_');
			sb.append(split[2]);
			sb.append('_');
			if(split[2].equals("sp.")){
				sb.append(split[3]);
				sb.append('_');
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			System.err.println(path);
			e.printStackTrace();
			return;
		}
		if(sb.length()>0){
			String name=f.getName();
			sb.append(name);
			f.renameTo(new File(sb.toString()));
		}
	}
	
	/*--------------------------------------------------------------*/

	private ArrayList<String> list=new ArrayList<String>();
	private PrintStream outstream=System.err;
	private static boolean verbose=false;
	
}
