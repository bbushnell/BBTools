package scalar;

import java.util.ArrayList;

import clade.Clade;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.LineParser1;
import sketch.Sketch;
import stream.Read;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntList;
import tracker.KmerTracker;

public class ScalarData{

	ScalarData(boolean storeTIDs, boolean storeNames, long numericID_) {
		if(storeTIDs) {taxIDs=new IntList();}
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
		while((line=bf.nextLine())!=null){
			if(line.length>0 && line[0]!='#'){
				bytesProcessed+=(line.length+1);
				pointsProcessed++;
				parser.set(line);
				gc.add(parser.parseFloat(0));
				hh.add(parser.parseFloat(1));
				caga.add(parser.parseFloat(2));
				if(parser.terms()>3 && taxIDs!=null) {taxIDs.add(parser.parseInt(3));}
				if(parser.terms()>4 && names!=null) {names.add(parser.parseString(4));}
			}
		}
		bf.close();
		return this;
	}
	
	public final void print(String outname,
			boolean header, boolean printTID, boolean printName) {
		FileFormat ffout=FileFormat.testOutput(outname, FileFormat.TXT, null, true, true, false, false);
		print(ffout, header, printTID, printName);
	}
	
	public final void print(FileFormat ffout,
			boolean header, boolean printTID, boolean printName) {
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);
		print(bsw, header, printTID, printName);
		bsw.poison();
	}
	
	public final void print(ByteStreamWriter bsw,
			boolean header, boolean printTID, boolean printName) {
		if(bsw==null) {return;}
		if(header) {
			bsw.print("#");
			bsw.print("GC\tHH\tCAGA");
			if(printTID) {bsw.print("\tTaxID");}
			if(printName) {bsw.print("\tName");}
			bsw.print('\n');
		}
		ByteBuilder bb=new ByteBuilder();
		FloatList[] fls=data();
		for(int i=0, len=fls[0].size(); i<len; i++) {
			for(int j=0; j<fls.length; j++) {
				bb.appendt(fls[j].get(i), 7);
			}
			if(taxIDs!=null && printTID) {bb.appendt(taxIDs.get(i));}
			if(names!=null && printName) {bb.appendt(names.get(i));}
			bb.replaceLast('\n');
			bsw.print(bb);
			bb.clear();
		}
	}

	public void add(Read r, KmerTracker dimers, 
			int interval, int minlen, int tid, boolean breakOnContig) {
		if(r==null) {return;}
		final byte[] bases=r.bases;
		if(clade!=null) {clade.add(bases, null);}
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
		if(taxIDs!=null) {taxIDs.add(tid);}
		if(names!=null) {names.add(name);}
		dimers.resetCount();
		pointsProcessed++;
	}
	
	public void add(ScalarData sd) {
		gc.append(sd.gc);
		hh.append(sd.hh);
		caga.append(sd.caga);
		if(taxIDs!=null) {taxIDs.addAll(sd.taxIDs);}
		if(names!=null) {names.addAll(sd.names);}
		readsProcessed+=sd.readsProcessed;
		basesProcessed+=sd.basesProcessed;
		bytesProcessed+=sd.bytesProcessed;
		pointsProcessed+=sd.pointsProcessed;
//		if(clade!=null) {clade.add(sd.clade);}
	}

	public int tid(int i){return taxIDs==null ? -1 : taxIDs.get(i);}
	public String name(int i){return names==null ? null : names.get(i);}

	public FloatList gc=new FloatList();
	public FloatList hh=new FloatList();
	public FloatList caga=new FloatList();
	public IntList taxIDs;
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

	public static boolean verbose=false;
	
}
