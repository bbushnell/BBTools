package prok;

import java.util.ArrayList;

import json.JsonObject;
import shared.Tools;
import structures.ByteBuilder;

public class ScoreTracker {
	
	public ScoreTracker(int type_){
		type=type_;
	}
	
	public void add(ScoreTracker st){
		geneStartScoreSum+=st.geneStartScoreSum;
		geneStopScoreSum+=st.geneStopScoreSum;
		geneInnerScoreSum+=st.geneInnerScoreSum;
		lengthSum+=st.lengthSum;
		
		geneStartScoreCount+=st.geneStartScoreCount;
		geneStopScoreCount+=st.geneStopScoreCount;
		geneInnerScoreCount+=st.geneInnerScoreCount;
		lengthCount+=st.lengthCount;
	}
	
	public void add(ArrayList<Orf>[] array){
		for(ArrayList<Orf> list : array){add(list);}
	}
	
	public void add(ArrayList<Orf> list){
		if(list==null){return;}
		for(Orf orf : list){
			if(orf.type==type){add(orf);}
		}
	}
	
	public void add(Orf orf){
		if(orf==null || orf.type!=type){return;}
		geneStartScoreSum+=orf.startScore;
		geneStopScoreSum+=orf.stopScore;
		geneInnerScoreSum+=orf.averageKmerScore();
		lengthSum+=orf.length();
		
		geneStartScoreCount++;
		geneStopScoreCount++;
		geneInnerScoreCount++;
		lengthCount++;
	}
	
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		bb.append("Start Score:          \t ").append(geneStartScoreSum/geneStartScoreCount, 4).nl();
		bb.append("Stop Score:           \t ").append(geneStopScoreSum/geneStopScoreCount, 4).nl();
		bb.append("Inner Score:          \t ").append(geneInnerScoreSum/geneInnerScoreCount, 4).nl();
		bb.append("Length:               \t ").append(lengthSum/(double)lengthCount, 2).nl();
		if(genomeSize>0){
			bb.append("Approx Genic Fraction:\t ").append(Tools.min(1.0, lengthSum/(double)genomeSize), 4).nl();
		}
		return bb.toString();
	}
	
	public JsonObject toJson(){
		JsonObject jo=new JsonObject();
		jo.addLiteral("Start Score", geneStartScoreSum/geneStartScoreCount, 4);
		jo.addLiteral("Stop Score", geneStopScoreSum/geneStopScoreCount, 4);
		jo.addLiteral("Inner Score", geneInnerScoreSum/geneInnerScoreCount, 4);
		jo.addLiteral("Length", lengthSum/(double)lengthCount, 2);
		if(genomeSize>0){
			jo.addLiteral("Approx Genic Fraction", Tools.min(1.0, lengthSum/(double)genomeSize), 4);
		}
		return jo;
	}
	
	long geneStartScoreCount=0;
	long geneStopScoreCount=0;
	long geneInnerScoreCount=0;
	long lengthCount=0;
	
	double geneStartScoreSum=0;
	double geneStopScoreSum=0;
	double geneInnerScoreSum=0;
	long lengthSum=0;
	
	long genomeSize=0;
	
	final int type;
	
}
