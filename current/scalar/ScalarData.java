package scalar;

import java.util.ArrayList;

import clade.Clade;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.LineParser1;
import sketch.Sketch;
import sketch.SketchMakerMini;
import stream.Read;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntList;
import tracker.KmerTracker;

public class ScalarData implements Comparable<ScalarData>{

	ScalarData(boolean storeNames, long numericID_){
		if(storeNames) {names=new ArrayList<String>();}
		numericID=numericID_;
	}
	
	public FloatList[] data() {return new FloatList[] {gc, hh, caga};}
	public FloatList[] reorder(int[] order) {
		FloatList[] list=new FloatList[] {gc, hh, caga};
		if(order!=null) {
			list=new FloatList[] {list[order[0]], list[order[1]], list[order[2]]};
		}
		return list;
	}
	
	public ScalarData readTSV(FileFormat ffin1){

		ByteFile bf=ByteFile.makeByteFile(ffin1);
		LineParser1 parser=new LineParser1('\t');

		byte[] line;
		String prev=null;
		while((line=bf.nextLine())!=null){
			if(line.length>0 && line[0]!='#'){
				bytesProcessed+=(line.length+1);
				pointsProcessed++;
				parser.set(line);
				gc.add(parser.parseFloat(0));
				hh.add(parser.parseFloat(1));
				caga.add(parser.parseFloat(2));
				if(parser.terms()>3) {taxIDs.add(parser.parseInt(3));}
				if(parser.terms()>4 && names!=null) {
					String name=parser.parseString(4);
					names.add(name!=null && name.equals(prev) ? prev : name);
					prev=name;
				}
			}
		}
		bf.close();
		return this;
	}
	
	public final void print(String outname, boolean header,
			boolean printName, boolean printPos, int interval) {
		FileFormat ffout=FileFormat.testOutput(outname, FileFormat.TXT, null, true, true, false, false);
		print(ffout, header, printName, printPos, interval);
	}
	
	public final void print(FileFormat ffout, boolean header, 
			boolean printName, boolean printPos, int interval) {
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);
		print(bsw, header, printName, printPos, interval);
		bsw.poison();
	}
	
	public static String header(boolean sideHeader, boolean printName, boolean printPos) {
		ByteBuilder bb=new ByteBuilder();
		bb.append("#");
		if(sideHeader) {bb.appendt("Header");}
		bb.append("GC\tHH\tCAGA");
		if(true) {bb.append("\tTaxID");}
		if(printPos) {bb.append("\tPos");}
		if(printName) {bb.append("\tName");}
		return bb.nl().toString();
	}
	
	public ByteBuilder mean(boolean sideHeader, String name) {		
		ByteBuilder bb=new ByteBuilder();
		if(sideHeader) {bb.appendt("Mean");}
		bb.appendt(gc.mean(),5);
		bb.appendt(hh.mean(),5);
		bb.appendt(caga.mean(),5);
		int tid=(taxIDs==null ? 0 : taxIDs.modeUnsorted());
		bb.appendt(tid);
		if(name!=null) {bb.appendt(name);}
		bb.replaceLast('\n');
		return bb;
	}
	
	public ByteBuilder stdev(boolean sideHeader, String name) {		
		ByteBuilder bb=new ByteBuilder();
		if(sideHeader) {bb.appendt("STDev");}
		bb.appendt(gc.stdev(),5);
		bb.appendt(hh.stdev(),5);
		bb.appendt(caga.stdev(),5);
		int tid=(taxIDs==null ? 0 : taxIDs.modeUnsorted());
		bb.appendt(tid);
		if(name!=null) {bb.appendt(name);}
		bb.replaceLast('\n');
		return bb;
	}
	
	public final void print(ByteStreamWriter bsw, boolean header, 
			boolean printName, boolean printPos, int interval) {
		if(bsw==null) {return;}
		if(header) {bsw.print(header(false, printName, printPos));}
		ByteBuilder bb=new ByteBuilder();
		FloatList[] fls=data();
		String prevName=null;
		int pos=0;
		for(int i=0, len=fls[0].size(); i<len; i++) {
			for(int j=0; j<fls.length; j++) {bb.appendt(fls[j].get(i), decimals);}
			if(taxIDs!=null){bb.append(taxIDs.get(i));}
			bb.tab();
			if(names!=null && printName) {
				String name=names.get(i);
				boolean match=(name!=null && name.equals(prevName));
				prevName=name;
				pos=match ? pos+interval : 0;
				if(printPos) {bb.appendt(pos);}
				if(printName) {bb.appendt(names.get(i));}
			}
			bb.replaceLast('\n');
			bsw.print(bb);
			bb.clear();
		}
	}

	public void add(Read r, KmerTracker dimers, SketchMakerMini smm,
			int interval, int minlen, int tid, boolean breakOnContig) {
		if(r==null) {return;}
		final byte[] bases=r.bases;
		if(makeClade && r.length()>=minCladeSize) {
			if(clade==null) {clade=new Clade(-1, -1, null);}
			clade.add(bases, null);
		}
		if(makeSketch && r.length()>=minSketchSize) {smm.processRead(r);}
		if(breakOnContig && r.length()<minlen) {return;}
		readsProcessed++;
		basesProcessed+=r.length();
		if(breakOnContig) {dimers.clearAll();}
		if(parseTID && tid<0) {tid=bin.BinObject.parseTaxID(r.name());}
		if(dimers.window>0) {
			for(byte b : bases) {
				boolean newValid=dimers.addWindowed(b);
				if(newValid && interval>0 && dimers.count()>=interval) {
					toInterval(dimers, tid, r.name());
				}
			}
		}else {dimers.add(bases);}
		if(verbose) {System.err.println("dimers.count()="+dimers.count()+", minlen="+minlen);}
		if(dimers.count()>=minlen) {toInterval(dimers, tid, r.name());}
	}

	private void toInterval(KmerTracker dimers, int tid, String name) {
		if(verbose) {System.err.println("calling toInterval(dimers)");}
		gc.add(dimers.GC());
		hh.add(dimers.HH());
		caga.add(dimers.CAGA());
		taxIDs.add(tid);
		if(names!=null) {names.add(name);}
		dimers.resetCount();
		pointsProcessed++;
	}
	
	public void add(ScalarData sd) {
		assert(sd!=this);
		gc.append(sd.gc);
		hh.append(sd.hh);
		caga.append(sd.caga);
		taxIDs.addAll(sd.taxIDs);
		if(names!=null) {names.addAll(sd.names);}
		readsProcessed+=sd.readsProcessed;
		basesProcessed+=sd.basesProcessed;
		bytesProcessed+=sd.bytesProcessed;
		pointsProcessed+=sd.pointsProcessed;
//		if(clade!=null) {clade.add(sd.clade);}
	}

	public int tid(int i){return taxIDs.get(i);}
	public String name(int i){return names==null ? null : names.get(i);}

	@Override
	public int compareTo(ScalarData o){
		return numericID<o.numericID ? -1 : numericID>o.numericID ? 1 : o.gc.size-gc.size;
	}

	public FloatList gc=new FloatList();
	public FloatList hh=new FloatList();
	public FloatList caga=new FloatList();
	public IntList taxIDs=new IntList();
	public ArrayList<String> names;

	long readsProcessed=0;
	long basesProcessed=0;
	long bytesProcessed=0;
	long pointsProcessed=0;
	
	Clade clade;
	Sketch sketch;
	
	final long numericID;
	
	public static boolean parseTID=false;
	public static boolean makeClade=false;
	public static boolean makeSketch=false;
	public static int minCladeSize=2000;
	public static int minSketchSize=5000;
	public static int decimals=5;

	public static boolean verbose=false;
	
}
