package prok;

import java.util.Arrays;

import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.LongHashSet;

class StatsContainer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	StatsContainer(int type_, int kInner, int framesInner, int kStart, int framesStart, int offsetStart, int kStop, int framesStop, int offsetStop){
		type=type_;
		name=ProkObject.typeStrings[type];
		setInner(kInner, framesInner);
		setStart(kStart, framesStart, offsetStart);
		setStop(kStop, framesStop, offsetStop);
	}

	StatsContainer(int type_){
		type=type_;
		name=ProkObject.typeStrings[type];
	}
	
	void setInner(int kInner, int framesInner){
		assert(inner==null);
		statsArray[0]=inner=new FrameStats(name+" inner", kInner, framesInner, 0);
	}
	
	void setStart(int kStart, int framesStart, int offsetStart){
		assert(start==null);
		statsArray[1]=start=new FrameStats(name+" start", kStart, framesStart, offsetStart);
	}
	
	void setStop(int kStop, int framesStop, int offsetStop){
		assert(stop==null);
		statsArray[2]=stop=new FrameStats(name+" stop", kStop, framesStop, offsetStop);
	}
	
	void setInner(FrameStats fs){
		assert(inner==null);
		assert(fs!=null);
		statsArray[0]=inner=fs;
	}
	
	void setStart(FrameStats fs){
		assert(start==null);
		statsArray[1]=start=fs;
	}
	
	void setStop(FrameStats fs){
		assert(stop==null);
		statsArray[2]=stop=fs;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public String toString(){
		return appendTo(new ByteBuilder()).toString();
	}
	
	public ByteBuilder appendTo(ByteBuilder bb) {
		bb.append("#name\t").append(name).nl();
		bb.append("#type\t").append(type).nl();
		bb.append("#count\t").append(lengthCount).nl();
		bb.append("#lengthSum\t").append(lengthSum).nl();
		//lengths
		bb.append("#contains\t").append(3).nl();
		if(ProkObject.processType(type)) {
			for(FrameStats fs : statsArray){
				fs.appendTo(bb);
			}
		}else{
			for(FrameStats fs : statsArray){
				fs.append0(bb);
			}
		}
		return bb;
	}
	
	public void clear(){
		for(int i=0; i<statsArray.length; i++){
			if(statsArray[i]!=null){
				statsArray[i].clear();
			}
		}

		assert(inner==statsArray[0]);
		assert(start==statsArray[1]);
		assert(stop==statsArray[2]);
		
		Arrays.fill(lengths, 0);
		lengthSum=0;
		lengthCount=0;
		calculate();
	}
	
	public void setFrom(StatsContainer sc){
		assert(sc.name.equals(name));
		for(int i=0; i<statsArray.length; i++){
			FrameStats fs=sc.statsArray[i];
			if(statsArray[i]==null){
				statsArray[i]=new FrameStats(fs.name, fs.k, fs.frames, fs.leftOffset);
				statsArray[i].add(fs);
			}else{
				statsArray[i].setFrom(fs);
			}
		}
		inner=statsArray[0];
		start=statsArray[1];
		stop=statsArray[2];
		
		for(int i=0; i<lengths.length; i++){lengths[i]=sc.lengths[i];}
		lengthSum=sc.lengthSum;
		lengthCount=sc.lengthCount;
		calculate();
	}
	
	public void add(StatsContainer sc){
		assert(sc.name.equals(name));
		for(int i=0; i<statsArray.length; i++){
			FrameStats fs=sc.statsArray[i];
			if(statsArray[i]==null){
				statsArray[i]=new FrameStats(fs.name, fs.k, fs.frames, fs.leftOffset);
			}
			statsArray[i].add(fs);
		}

		inner=statsArray[0];
		start=statsArray[1];
		stop=statsArray[2];
		
		Tools.add(lengths, sc.lengths);
		lengthSum+=sc.lengthSum;
		lengthCount+=sc.lengthCount;
		calculate();
	}
	
	public void multiplyBy(double mult) {
		for(int i=0; i<statsArray.length; i++){
			FrameStats fs=statsArray[i];
			fs.multiplyBy(mult);
		}
		
		Tools.multiplyBy(lengths, mult);
		lengthSum=Math.round(lengthSum*mult);
		lengthCount=Math.round(lengthCount*mult);
		calculate();
	}
	
	public void calculate(){
		for(int i=0; i<statsArray.length; i++){
			statsArray[i].calculate();
		}
		lengthAvg=(int)(lengthSum/Tools.max(1.0, lengthCount));
		invLengthAvg=1f/Tools.max(1, lengthAvg);
	}
	
	public void addLength(int x){
		lengthSum+=x;
		lengthCount++;
		lengths[Tools.min(x, lengths.length-1)]++;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	FrameStats inner;
	FrameStats start;
	FrameStats stop;
	final FrameStats[] statsArray=new FrameStats[3];
	
	int kLongLen(){return ProkObject.kLongLen(type);}
	LongHashSet kmerSet(){return ProkObject.kmerSet(type);}
	Read[] consensusSequence(){return ProkObject.consensusReads(type);}
	float minIdentity(){return ProkObject.minID(type);}

	public int startSlop(){return ProkObject.startSlop(type);}
	public int stopSlop(){return ProkObject.stopSlop(type);}
	
	final String name;
	long lengthSum=0;
	long lengthCount=0;
	int lengthAvg=-1;
	float invLengthAvg;
	
	int[] lengths=new int[5000];
	public final int type;
	
}
