package kmer;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import shared.KillSwitch;
import shared.Shared;
import shared.Tools;

public class OwnershipThread extends Thread {
	
	public static void clear(AbstractKmerTable[] tables){
		process(tables, CLEAR);
	}
	
	public static void initialize(AbstractKmerTable[] tables){
		process(tables, INITIALIZE);
	}
	
	private static void process(AbstractKmerTable[] tables, int mode){
		if(tables.length<2){
			if(mode==INITIALIZE){
				for(AbstractKmerTable akt : tables){akt.initializeOwnership();}
			}else if(mode==CLEAR){
				for(AbstractKmerTable akt : tables){akt.clearOwnership();}
			}else{
				KillSwitch.kill("Bad mode: "+mode);
			}
			return;
		}
		final int threads=Tools.min(Shared.threads(), tables.length);
		final AtomicInteger next=new AtomicInteger(0);
		ArrayList<OwnershipThread> alpt=new ArrayList<OwnershipThread>(threads);
		for(int i=0; i<threads; i++){alpt.add(new OwnershipThread(tables, mode, next));}
		for(OwnershipThread pt : alpt){pt.start();}
		
		for(OwnershipThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	public OwnershipThread(AbstractKmerTable[] tables_, int mode_, AtomicInteger next_){
		tables=tables_;
		mode=mode_;
		next=next_;
	}
	
	@Override
	public void run(){
		for(int i=next.getAndIncrement(); i<tables.length; i=next.getAndIncrement()){
			if(mode==INITIALIZE){
				tables[i].initializeOwnership();
			}else if(mode==CLEAR){
				tables[i].clearOwnership();
			}else{
				KillSwitch.kill("Bad mode: "+mode);
			}
		}
	}
	
	private final AbstractKmerTable[] tables;
	private final AtomicInteger next;
	private final int mode;

	public static final int INITIALIZE=0;
	public static final int CLEAR=1;
	
}
