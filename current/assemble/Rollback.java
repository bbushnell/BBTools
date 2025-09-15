package assemble;

import stream.Read;
import structures.IntList;

public class Rollback {

	public Rollback(Read r){
		this(r, null);
	}

	public Rollback(Read r, IntList counts){
		id0=r.id;
		flags0=r.flags;
		bases0=r.bases.clone();
		quals0=(r.quality==null ? null : r.quality.clone());
		counts0=(counts==null ? null : counts.copy());
	}
	
	public void rollback(Read r){
		rollback(r, null);
	}
	
	public void rollback(Read r, IntList counts){
		r.id=id0;
		r.flags=flags0;
		if(r.length()==bases0.length){
			System.arraycopy(bases0, 0, r.bases, 0, bases0.length);
			if(quals0!=null){System.arraycopy(quals0, 0, r.quality, 0, quals0.length);}
			if(counts!=null){System.arraycopy(counts0.array, 0, counts.array, 0, counts0.size);}
		}else{
			r.bases=bases0;
			r.quality=quals0;
			if(counts!=null){
				counts.clear();
				counts.addAll(counts0);
			}
		}
	}
	
	final String id0;
	final int flags0;
	final byte[] bases0, quals0;
	public final IntList counts0;
	
}
