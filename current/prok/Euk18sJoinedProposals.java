package prok;

import java.util.ArrayList;
import fileIO.ByteFile;
import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;

/** Frozen proposal replay for the D45 fast-caller experiment. Loads only geometry
 * and model identity; saved truth-containment and role fields never affect calls.
 * The caller performs MSA/identity acceptance and its normal final path selection.
 * @author Brian Bushnell, Raiden */
final class Euk18sJoinedProposals {
	static Euk18sJoinedProposals load(String path,String[] names,int support){
		require(path!=null && names!=null && (support==1 || support==2),"Joined replay requires an explicit table, actual model names and side support1 or2");
		final Euk18sJoinedProposals result=new Euk18sJoinedProposals();final ObjectIntMap<String> models=new ObjectIntMap<String>(String.class);
		for(int i=0;i<names.length;i++){require(models.put(names[i],i)==-1,"Unique model names bind proposal alignments");}
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,HEADER,path);
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==19,"Preserve the phase1 proposal schema");
				final int arm=p.parseInt(3);require(arm==1 || arm==2,"Unexpected phase1 side-support arm");if(arm!=support){continue;}
				for(int term:new int[]{14,16,17}){require(p.termEquals("true",term) || p.termEquals("false",term),"Proposal eligibility/tandem/reset must be explicit booleans");}
				if(!p.termEquals("true",14) || p.termEquals("true",16) || p.termEquals("true",17)){continue;}
				final String target=p.parseString(0),strand=p.parseString(2),model=p.parseString(4);
				require((strand.equals("+") || strand.equals("-")) && models.contains(model),"Proposal orientation and model must bind an actual loaded consensus");
				final int start=p.parseInt(6),stop=p.parseInt(7);require(start>0 && stop>=start,"Joined proposals retain physical one-based inclusive bounds");
				final String key=target+"\t"+strand,unique=key+"\t"+model+"\t"+start+"\t"+stop;
				if(!seen.add(unique)){continue;}// Same window/model at another site is the same MSA.
				Group group=result.groups.get(key);if(group==null){group=new Group();result.groups.put(key,group);}
				group.rows.add(new Proposal(start,stop,models.get(model)));result.windows++;
			}
		}finally{require(!in.close(),"Missing proposal rows cannot silently change the replay comparison");}
		System.err.println("EUK18S_JOINED_PROPOSALS support="+support+" distinctWindowModels="+result.windows+" savedGeometryOnly=true");return result;
	}
	Group get(String name,int strand){require(strand==0 || strand==1,"Caller strand must be0or1");return groups.get(RrnaResourceIO.first(name)+"\t"+(strand==0?"+":"-"));}
	static final class Proposal {
		Proposal(int start_,int stop_,int model_){assert(start_>0 && stop_>=start_ && model_>=0):"A joined window must retain its physical interval and model";start=start_;stop=stop_;model=model_;}
		final int start,stop,model;
	}
	static final class Group {final ArrayList<Proposal> rows=new ArrayList<Proposal>();}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	final ObjectMap<String,Group> groups=new ObjectMap<String,Group>(String.class,Group.class);int windows;
	static final String HEADER="shred\trole\tstrand\tsideSupport\tmodel\tRFboundary\tstart1\tstop1\torientedRawStart0\torientedRawStop0\tgapLow\tgapHigh\tleftKeys\trightKeys\teligible\tambiguous\ttandem\tmodelReset\tcontainsPanelTruth";
}
