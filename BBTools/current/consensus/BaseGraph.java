package consensus;

import java.io.Serializable;
import java.util.Arrays;

import aligner.Aligner;
import aligner.AlignmentResult;
import aligner.FlatAligner2;
import dna.AminoAcid;
import prok.GeneCaller;
import shared.KillSwitch;
import shared.Tools;
import shared.TrimRead;
import stream.FASTQ;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;

/**
 * A graph of the possible alignments of a reference sequence.
 * 
 * @author Brian Bushnell
 * @date September 6, 2019
 *
 */
public class BaseGraph extends ConsensusObject implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 3000198119165515066L;
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public BaseGraph(String name_, byte[] bases_, byte[] quals_, long numericID_, int pad_){
		name=name_;
		pad=pad_;
		original=pad(bases_, pad, (byte)'N');
		refWeights=makeWeightsFromQual(quals_);
		numericID=numericID_;

		ref=new BaseNode[original.length];
		del=new BaseNode[original.length];
		for(int i=0; i<original.length; i++){
			byte b=original[i];
			ref[i]=new BaseNode(b, REF, i);
			del[i]=new BaseNode(b, DEL, i);
		}
		for(int i=0; i<ref.length-1; i++){
			ref[i].refEdge=ref[i+1];
		}
		for(int i=0; i<ref.length; i++){
			ref[i].add(original[i], 1);
		}
	}
	
	private float[] makeWeightsFromQual(byte[] quals){
		float[] w=new float[original.length];
		Arrays.fill(w, 1f);
		for(int i=0; i<pad; i++){
			w[i]=w[w.length-i-1]=0;
		}
		if(quals!=null){
			float mult=1f/39;
			for(int i=0; i<quals.length; i++){
				int q=(Tools.max(0, quals[i]-2));
				w[i+pad]=0.6f+(0.4f*mult*q);
			}
		}
		return w;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	private static byte[] pad(byte[] in, int pad, byte symbol){
		if(pad<1){return in;}
		byte[] out=new byte[in.length+2*pad];
		Arrays.fill(out, symbol);
		for(int i=0; i<in.length; i++){
			out[i+pad]=in[i];
		}
		return out;
	}
	
	/** This method should be threadsafe */
	public void add(Read r){
		add(r, 0, 0);
	}
	
	/** This method should be threadsafe */
	public void add(Read r, int leftBorder, int rightBorder){
		final int maxQpos=r.length()-rightBorder;
		r.toLongMatchString(false);
		final byte[] match=r.match;
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		assert(match!=null && bases!=null) : r;
		final int start=r.start;
		int qpos=0, rpos=start;
		byte prevState='?';
		BaseNode prevNode=(rpos<=0 ? null : ref[rpos-1]);
		
		final int mapq=(r.samline==null ? 20 : r.samline.mapq+1);
		int msSum=0;
		int delSum=0;
		int insSum=0;
		long qualSum=0;
		int identity=(int)(r.identity()*100);
		for(int mpos=0; mpos<match.length && rpos<ref.length; mpos++){
			final byte m=match[mpos];
			final BaseNode next;
			final int q;
			assert(qpos<bases.length) : "\n"+r+"\n"+r.samline+"\n"+r.length()+", "+r.start+", "+r.stop+", "+
				r.samline.calcCigarLength(true, false)+", "+r.samline.calcCigarReadLength(true, false);
			final byte b=bases[qpos];
			
//			System.err.println("mpos="+mpos+", qpos="+qpos+", rpos="+rpos+", b="+Character.toString((char)b)+", m="+Character.toString((char)m)+", prev="+prevNode);
			
			if(m=='m' || m=='S' || m=='N'){
				next=(rpos<0 ? null : ref[rpos]);
				q=(quals==null ? fakeQuality : quals[qpos]);

				qpos++;
				rpos++;
				msSum++;
				qualSum+=q;
			}else if(m=='D'){
				next=(rpos<0 ? null : del[rpos]);
				q=(quals==null ? fakeQuality : (quals[qpos]+quals[qpos+1])/2);

				rpos++;
				delSum++;
			}else if(m=='I'){
				assert(prevNode!=null) : new String(r.match+"\n"+r.samline); //Alignments must not start with I.
				synchronized(prevNode) {
					if(prevNode.insEdge==null){
						prevNode.insEdge=new BaseNode('.', INS, rpos);
					}
					next=prevNode.insEdge;
				}
				q=(quals==null ? fakeQuality : quals[qpos]);

				qpos++;
				insSum++;
				qualSum+=q;
			}else if(m=='C'){
				next=null;
				
				q=(quals==null ? fakeQuality : quals[qpos]);
				qpos++;
				rpos++;
			}else if(m=='X' || m=='Y'){
				assert(rpos<0 || rpos>=ref.length) : rpos+", "+ref.length+", "+r.start+"\n"+new String(bases)+"\n"+new String(original)+"\n"+new String(match);
				next=null;
				
				q=(quals==null ? fakeQuality : quals[qpos]);
				qpos++;
				rpos++;
			}else{
				next=null;
				
				q=1;
				assert(false) : Character.toString((char)m)+"\n"+new String(match)+"\n"+r.id;
			}
			
			if(next!=null && qpos>=leftBorder && qpos<=maxQpos){
				assert(qpos>=leftBorder && qpos<=maxQpos):qpos+", "+leftBorder+", "+maxQpos+", "+r.length();
//				System.err.println("Adding "+Character.toString(b)+" to "+next);
				int weight=q+1;
				if(useMapq){
					weight=(int)(Math.ceil(Math.sqrt(q*mapq)));
				}
				if(invertIdentity){
					weight=weight+Tools.max(0, identityCeiling-identity);
				}
				synchronized(next){
					next.add(b, weight);
				}
			}
			prevNode=next;
			prevState=m;
		}
		
		synchronized(this){
			readTotal++;
			baseTotal+=bases.length;
			symbolTotal+=match.length;
			msTotal+=msSum;
			insTotal+=insSum;
			delTotal+=delSum;
			qualTotal+=qualSum;
		}
	}
	
	/** 
	 * Generates a score indicating how well the alignment matches the model.
	 * Scores range from 1 to -1 times the alignment length; a random sequence
	 * might be around 0.  Less-conserved models result in lower scores,
	 * even for optimally-matching alignments. 
	 * A relative score divides by the alignment length. */
	public float scoreOld(Read r, boolean local, boolean relative){
		float sum=0;
		float minSum=0;
		float maxSum=0;
		float peak=0;
		float minPeak=0;
		float maxPeak=0;
		
//		new Exception().printStackTrace();
		
		r.toLongMatchString(false);
		final byte[] match=r.match;
		final byte[] bases=r.bases;
		assert(match!=null && bases!=null);
		final int start=r.start;
		int qpos=0, rpos=start;
		byte prevState='?';
		BaseNode prevNode=(rpos<=0 ? null : ref[rpos-1]);
		
		int msSum=0;
		int delSum=0;
		int insSum=0;
		for(int mpos=0; mpos<match.length && rpos<ref.length; mpos++){
			final byte m=match[mpos];
			final BaseNode next;
			assert(qpos<bases.length) : "\n"+r+"\n"+r.samline+"\n"+r.length()+", "+r.start+", "+r.stop+", "+
				r.samline.calcCigarLength(true, false)+", "+r.samline.calcCigarReadLength(true, false);
			final byte b=bases[qpos];

//			System.err.println("mpos="+mpos+", qpos="+qpos+", rpos="+rpos+", b="+Character.toString((char)b)+", m="+Character.toString((char)m)+", prev="+prevNode);
//			System.err.println("mpos="+mpos+", qpos="+qpos+", rpos="+rpos+", b="+Character.toString((char)b)+", m="+Character.toString((char)m)+", prev="+prevNode);
			
			if(m=='m' || m=='S' || m=='N'){
				next=(rpos<0 || rpos>=ref.length ? null : ref[rpos]);
				if(verbose){
					next.delEdge=(rpos+1>=del.length ? null : del[rpos+1]);
					System.err.println("mpos="+mpos+", qpos="+qpos+", rpos="+rpos+", b="+Character.toString((char)b)+
							", m="+Character.toString((char)m)+", next="+next);
				}
				assert(next!=null) : "Ref events should never be null. rpos="+rpos+", start="+r.start+"\n"+new String(match);

				float baseScore=next.baseScore(b);
				float minBaseScore=next.baseScore(next.minBase());
				float maxBaseScore=next.baseScore(next.maxBase());
				float stateProb=next.countSum/(Tools.max(1f, countSum(rpos)));
				float score=baseScore*stateProb+0.1f*(stateProb-1f);
				float minScore=minBaseScore*stateProb+0.1f*(stateProb-1f);
				float maxScore=maxBaseScore*stateProb+0.1f*(stateProb-1f);
//				float score=(stateProb-0.5f)+0.5f*(baseScore);
				sum+=score;
				minSum+=minScore;
				maxSum+=maxScore;
				if(verbose){
					System.err.println("countSum="+next.countSum+", countSum(rpos)="+countSum(rpos)+", "
							+ "baseScore="+baseScore+", stateProb="+stateProb+"\n"
									+ "score="+score+", minScore="+minScore+", minScore="+maxScore+", sum="+sum+", minSum="+minSum+", maxSum="+maxSum);
					System.err.println();
				}
				
				qpos++;
				rpos++;
				msSum++;
			}else if(m=='D'){
				next=(rpos<0 || rpos>=ref.length ? null : del[rpos]);
				if(verbose){
					next.delEdge=(rpos+1>=del.length ? null : del[rpos+1]);
					next.refEdge=(rpos+1>=ref.length ? null : ref[rpos+1]);
					System.err.println("mpos="+mpos+", qpos="+qpos+", rpos="+rpos+", b="+Character.toString((char)b)+
							", m="+Character.toString((char)m)+", next="+next);
				}
				assert(next!=null) : "Del events should never be null. rpos="+rpos+", start="+r.start+"\n"+new String(match);
				float stateProb=next.countSum/(Tools.max(1f, countSum(rpos)));
				float score=2*(stateProb-0.5f);
//				assert(score<=0) : score;//This is fine if it comes from a non-consensus; del can be the dominant mode after all.
				sum+=score;
				minSum+=score;
				if(verbose){
					System.err.println("countSum="+next.countSum+", countSum(rpos)="+countSum(rpos)+", "
							+ "stateProb="+stateProb+", score="+score+", sum="+sum);
					System.err.println();
				}
				
				rpos++;
				delSum++;
			}else if(m=='I'){
				assert(prevNode==null || prevNode.type==REF || prevNode.type==INS) : prevNode;
				if(prevNode==null){next=null;}
				else if(prevNode.insEdge!=null){next=prevNode.insEdge;}
				else if(prevNode.type==INS){next=prevNode;}
				else{next=null;}
				if(verbose){
					if(next!=null){next.refEdge=(rpos+1>=ref.length ? null : ref[rpos+1]);}
					System.err.println("mpos="+mpos+", qpos="+qpos+", rpos="+rpos+", b="+Character.toString((char)b)+
							", m="+Character.toString((char)m)+", next="+next);
				}

				//TODO: Effect of base score should be modified by stateprob
				final float baseScore, stateProb;
				final float minBaseScore, maxBaseScore;
				if(next==null){
					baseScore=-1;
					minBaseScore=maxBaseScore=-1;
					stateProb=0;
				}else{
					baseScore=next.baseScore(b);
					minBaseScore=next.baseScore(next.minBase());
					maxBaseScore=next.baseScore(next.maxBase());
					stateProb=next.countSum/(Tools.max(1f, countSum(rpos)));
				}
				float score=(stateProb-0.5f)+(0.25f*(baseScore-1));
				float minScore=(stateProb-0.5f)+(0.25f*(minBaseScore-1));
				float maxScore=(stateProb-0.5f)+(0.25f*(maxBaseScore-1));
				sum+=score;
				minSum+=minScore;
				maxSum+=maxScore;
				if(verbose){
					System.err.println("countSum="+(next==null ? 0 : next.countSum)+", countSum(rpos)="+countSum(rpos)+", "
							+ "baseScore="+baseScore+", stateProb="+stateProb+", score="+score+", sum="+sum);
					System.err.println();
				}
				
				qpos++;
				insSum++;
			}else if(m=='C'){
				next=null;
				
				qpos++;
				rpos++;
			}else if(m=='X' || m=='Y'){
				assert(rpos<0 || rpos>=ref.length) : rpos+", "+ref.length+", "+r.start+"\n"+new String(bases)+"\n"+new String(original)+"\n"+new String(match);
				next=null;
				
				qpos++;
				rpos++;
			}else{
				next=null;
				
				assert(false) : Character.toString((char)m)+"\n"+new String(match)+"\n"+r.id;
			}
			
			prevNode=next;
			prevState=m;
			if(local){
				if(sum<0){
					sum=0;
					minSum=0;
					maxSum=0;
				}else if(sum>peak){
					peak=sum;
					minPeak=minSum;
					maxPeak=maxSum;
				}
			}
		}
		
		boolean oldMode=false;
		if(oldMode){
			float score=local ? peak : sum;
			int length=msSum+delSum+insSum;
			return relative ? (1+(score/length))*0.5f : score;//Relative score does not play well with local alignments
		}else{
//			System.err.println(minSum+", "+maxSum+", "+sum);
			float a=local ? peak : sum;
			float b=local ? minPeak : minSum;
			float c=local ? maxPeak : maxSum;
			return relative ? (a-b)/(c-b) : a; //This changes to a 0-1 scale of worst to best possible alignment
		}
	}
	
	/** 
	 * Generates a score indicating how well the alignment matches the model.
	 * Scores range from 1 to -1 times the alignment length; a random sequence
	 * might be around 0.  Less-conserved models result in lower scores,
	 * even for optimally-matching alignments. 
	 * A relative score divides by the alignment length. */
	public float score(Read r, boolean local, boolean relative){
		float sum=0;
		float minSum=0;
		float maxSum=0;
		float peak=0;
		float minPeak=0;
		float maxPeak=0;
		
//		new Exception().printStackTrace();
		
		r.toLongMatchString(false);
		final byte[] match=r.match;
		final byte[] bases=r.bases;
		assert(match!=null && bases!=null);
		final int start=r.start;
		int qpos=0, rpos=start;
		byte prevState='?';
		BaseNode prevNode=(rpos<=0 ? null : ref[rpos-1]);
		
		int msSum=0;
		int delSum=0;
		int insSum=0;
		for(int mpos=0; mpos<match.length && rpos<ref.length; mpos++){
			final byte m=match[mpos];
			final BaseNode next;
			assert(qpos<bases.length) : "\n"+r+"\n"+r.samline+"\n"+r.length()+", "+r.start+", "+r.stop+", "+
				r.samline.calcCigarLength(true, false)+", "+r.samline.calcCigarReadLength(true, false);
			final byte b=bases[qpos];

//			System.err.println("mpos="+mpos+", qpos="+qpos+", rpos="+rpos+", b="+Character.toString((char)b)+", m="+Character.toString((char)m)+", prev="+prevNode);
//			System.err.println("mpos="+mpos+", qpos="+qpos+", rpos="+rpos+", b="+Character.toString((char)b)+", m="+Character.toString((char)m)+", prev="+prevNode);
			
			if(m=='m' || m=='S' || m=='N'){
				next=(rpos<0 || rpos>=ref.length ? null : ref[rpos]);
				if(verbose){
					next.delEdge=(rpos+1>=del.length ? null : del[rpos+1]);
					System.err.println("mpos="+mpos+", qpos="+qpos+", rpos="+rpos+", b="+Character.toString((char)b)+
							", m="+Character.toString((char)m)+", next="+next);
				}
				assert(next!=null) : "Ref events should never be null. rpos="+rpos+", start="+r.start+"\n"+new String(match);

				float baseScore=next.baseScore(b);
				float minBaseScore=next.baseScore(next.minBase());
				float maxBaseScore=next.baseScore(next.maxBase());
				float stateProb=next.countSum/(Tools.max(1f, countSum(rpos)));
				float score=baseScore*stateProb+0.1f*(stateProb-1f);
				float minScore=minBaseScore*stateProb+0.1f*(stateProb-1f);
				float maxScore=maxBaseScore*stateProb+0.1f*(stateProb-1f);
//				float score=(stateProb-0.5f)+0.5f*(baseScore);
				sum+=score;
				minSum+=minScore;
				maxSum+=maxScore;
				if(verbose){
					System.err.println("countSum="+next.countSum+", countSum(rpos)="+countSum(rpos)+", "
							+ "baseScore="+baseScore+", stateProb="+stateProb+"\n"
									+ "score="+score+", minScore="+minScore+", minScore="+maxScore+", sum="+sum+", minSum="+minSum+", maxSum="+maxSum);
					System.err.println();
				}
				
				qpos++;
				rpos++;
				msSum++;
			}else if(m=='D'){
				next=(rpos<0 || rpos>=ref.length ? null : del[rpos]);
				if(verbose){
					next.delEdge=(rpos+1>=del.length ? null : del[rpos+1]);
					next.refEdge=(rpos+1>=ref.length ? null : ref[rpos+1]);
					System.err.println("mpos="+mpos+", qpos="+qpos+", rpos="+rpos+", b="+Character.toString((char)b)+
							", m="+Character.toString((char)m)+", next="+next);
				}
				assert(next!=null) : "Del events should never be null. rpos="+rpos+", start="+r.start+"\n"+new String(match);
				float stateProb=next.countSum/(Tools.max(1f, countSum(rpos)));
				float score=2*(stateProb-0.5f);
//				score=Tools.mid(-1, (score-0.1f)*1.1f, 0);
				float maxScore=Tools.max(0, score);
//				assert(score<=0) : score;//This is fine if it comes from a non-consensus; del can be the dominant mode after all.
				sum+=score;
				minSum+=score;
				maxSum+=maxScore;
				if(verbose){
					System.err.println("countSum="+next.countSum+", countSum(rpos)="+countSum(rpos)+", "
							+ "stateProb="+stateProb+", score="+score+", sum="+sum);
					System.err.println();
				}
				
				rpos++;
				delSum++;
			}else if(m=='I'){
				assert(prevNode==null || prevNode.type==REF || prevNode.type==INS) : prevNode;
				if(prevNode==null){next=null;}
				else if(prevNode.insEdge!=null){next=prevNode.insEdge;}
				else if(prevNode.type==INS){next=prevNode;}
				else{next=null;}
				if(verbose){
					if(next!=null){next.refEdge=(rpos+1>=ref.length ? null : ref[rpos+1]);}
					System.err.println("mpos="+mpos+", qpos="+qpos+", rpos="+rpos+", b="+Character.toString((char)b)+
							", m="+Character.toString((char)m)+", next="+next);
				}

				//TODO: Effect of base score should be modified by stateprob
				final float baseScore, stateProb;
				final float minBaseScore, maxBaseScore;
				if(next==null){
					baseScore=-1;
					minBaseScore=maxBaseScore=-1;
					stateProb=0;
				}else{
					baseScore=next.baseScore(b);
					minBaseScore=next.baseScore(next.minBase());
					maxBaseScore=next.baseScore(next.maxBase());
					stateProb=next.countSum/(Tools.max(1f, countSum(rpos)));
				}
				float score=(stateProb-0.5f)+(0.25f*(baseScore-1));
				float minScore=(stateProb-0.5f)+(0.25f*(minBaseScore-1));
				float maxScore=(stateProb-0.5f)+(0.25f*(maxBaseScore-1));
				
//				score=Tools.mid(-1, (score-0.1f)*1.1f, 0);
				minScore=Tools.min(score, minScore);
				maxScore=Tools.max(0, score, maxScore);
				
				sum+=score;
				minSum+=minScore;
				maxSum+=maxScore;
				if(verbose){
					System.err.println("countSum="+(next==null ? 0 : next.countSum)+", countSum(rpos)="+countSum(rpos)+", "
							+ "baseScore="+baseScore+", stateProb="+stateProb+", score="+score+", sum="+sum);
					System.err.println();
				}
				
				qpos++;
				insSum++;
			}else if(m=='C'){
				next=null;
				
				qpos++;
				rpos++;
			}else if(m=='X' || m=='Y'){
				assert(rpos<0 || rpos>=ref.length) : rpos+", "+ref.length+", "+r.start+"\n"+new String(bases)+"\n"+new String(original)+"\n"+new String(match);
				next=null;
				
				qpos++;
				rpos++;
			}else{
				next=null;
				
				assert(false) : Character.toString((char)m)+"\n"+new String(match)+"\n"+r.id;
			}
			
			prevNode=next;
			prevState=m;
			if(local){
				if(sum<0){
					sum=0;
					minSum=0;
					maxSum=0;
				}else if(sum>peak){
					peak=sum;
					minPeak=minSum;
					maxPeak=maxSum;
				}
			}
		}
		
		boolean oldMode=false;
		if(oldMode){
			float score=local ? peak : sum;
			int length=msSum+delSum+insSum;
			return relative ? (1+(score/length))*0.5f : score;//Relative score does not play well with local alignments
		}else{
//			System.err.println(minSum+", "+maxSum+", "+sum);
			float a=local ? peak : sum;
			float b=local ? minPeak : minSum;
			float c=local ? maxPeak : maxSum;
			return relative ? (a-b)/(c-b) : a; //This changes to a 0-1 scale of worst to best possible alignment
		}
	}
	
	/** Count of alignments over this reference position. */
	int countSum(int rpos){
		assert(rpos>=0 && rpos<ref.length);
		if(rpos<0 || rpos>=ref.length){return 0;}
		BaseNode dnode=del[rpos];
		BaseNode rnode=ref[rpos];
		int sum=(dnode==null ? 0 : dnode.countSum)+(rnode==null ? 0 : rnode.countSum);
		return sum;
	}
	
	@Override
	public ByteBuilder toText() {
		System.err.println("...");
		ByteBuilder bb=new ByteBuilder();
		bb.append(name).append(':');
		bb.append('{');
		if(ref.length>0){bb.append(ref[0].toText());}
		for(int i=0; i<ref.length; i++){
			if(i>0){bb.comma();}
			bb.append(ref[i].toText());
		}
		bb.append('}');
		return bb;
	}
	
	int maxDepth(){
		int maxDepth=0;
		for(int i=0; i<ref.length; i++) {
			final BaseNode dnode=del[i];
			final BaseNode rnode=ref[i];
			final int dc=dnode.countSum;
			final int rc=rnode.countSum;
			final int depth=dc+rc;
			maxDepth=Tools.max(maxDepth, depth);
		}
		return maxDepth;
	}
	
	public Read traverse() {
		//TODO: Note, this mutates the object due to refCount, subCount, insCount, delCount
		//Maybe they should be cleared first, or just removed if not needed.
		final ByteBuilder bb=new ByteBuilder();
		final ByteBuilder bq=new ByteBuilder();
		final IntList depthList=new IntList();
		
		int maxDepth=0;	
		
		subCount=0;
		refCount=0;
		delCount=0;
		insCount=0;
		
//		System.err.println("rpos\tdw\trw\tiw");
		final byte[] temp=KillSwitch.allocByte1D(2);
		for(int i=0; i<ref.length; i++) {
			
			final BaseNode dnode=del[i];
			final BaseNode rnode=ref[i];
			BaseNode inode=(rnode.insEdge==null ? dummy : rnode.insEdge);
			
			final int dw=dnode.weightSum;
			final int rw=rnode.weightSum;

			final int dc=dnode.countSum;
			final int rc=rnode.countSum;
			final int depth=dc+rc;
			maxDepth=Tools.max(maxDepth, depth);
			
			final float afMult=1f/(dc+rc);
			final float daf=dc*afMult;
			
			long weightSum=dw+rw;
			
			
//			System.err.println(i+"\t"+dw+"\t"+rw+"\t"+iw);
			
			if(rw>=dw || daf<MAF_del || noIndels){//Common case
				{
					rnode.consensus(temp);
					byte b=temp[0];
					byte q=temp[1];
					bb.append(b);
					depthList.add(depth);
//					System.err.println(i+": "+Character.toString(b));
					assert(b!='.') : rnode+", "+rnode.weightSum+", "+weightSum;
					if(qualityByMS){
						float prob=rnode.alleleDif();
						byte quality=(byte)Tools.mid(2, 41, Math.round(39*prob));
						bq.append(quality);
//						System.err.println("quality="+quality+" for "+Arrays.toString(rnode.acgtCount));
					}else{
						double quality=Tools.mid(2, 41, 10*Math.log10(rw/Tools.max(0.01f, weightSum-rw)));
						bq.append((byte)(Tools.min(q, quality)+FASTQ.ASCII_OFFSET));
					}
					if(b==rnode.refBase){
						refCount++;
					}else{
						subCount++;
					}
				}
				
				//Then add inodes
				while(inode!=null && !noIndels && inode.weightSum>=(weightSum-inode.weightSum) && inode.countSum*afMult>=MAF_ins){
					inode.consensus(temp);
					byte b=temp[0];
					byte q=temp[1];
					assert(b!='.') : inode+", "+inode.weightSum+", "+weightSum;
					bb.append(b);
					depthList.add(depth);
//					System.err.println(i+": "+Character.toString(b));
					if(qualityByMS){
						byte quality=(byte)Tools.mid(2, 41, Math.round(39*inode.alleleDif()));
						bq.append(quality);
					}else{
						double quality=Tools.mid(2, 41, 10*Math.log10(inode.weightSum/(rw+dw)));
						bq.append((byte)(Tools.min(q, quality)+FASTQ.ASCII_OFFSET));
					}
					insCount++;
					inode=inode.insEdge;
				}
			}else{
				//Deletion, do nothing
				delCount++;
			}
		}
		
		Read r=new Read(bb.toBytes(), bq.toBytes(), name, numericID);
		
		if(trimDepthFraction>0 || trimNs){
			final int trimDepth=Tools.max(1, (int)(trimDepthFraction*maxDepth));
			int left=0, right=0;
			while(left<depthList.size && (depthList.get(left)<trimDepth || r.bases[left]=='N')){left++;}
			while(right<depthList.size && (depthList.get(depthList.size-right-1)<trimDepth || r.bases[depthList.size-right-1]=='N')){right++;}
			if(left>0 || right>0){
				TrimRead.trimByAmount(r, left, right, 1);
			}
//			System.err.println(left+", "+right);
		}
		
		return r;
	}
	
	public void makeWeights() {
		
		refWeights=new float[ref.length];
//		insWeights=new float[ref.length];
//		delWeights=new float[ref.length];
		
		for(int i=0; i<ref.length; i++) {
			
			final BaseNode dnode=del[i];
			final BaseNode rnode=ref[i];
			BaseNode inode=(rnode.insEdge==null ? dummy : rnode.insEdge);
			
			final int dw=dnode.weightSum;
			final int rw=rnode.weightSum;

			final int dc=dnode.countSum;
			final int rc=rnode.countSum;
			final int depth=dc+rc;
			
			final float afMult=1f/depth;
			final float daf=dc*afMult;
			final float iaf=inode.countSum*afMult;
			
			long weightSum=dw+rw;
			

			refWeights[i]=0.5f*rnode.alleleDif()+0.5f;
//			insWeights[i]=0.5f*(1-iaf)+0.5f;
//			delWeights[i]=0.5f*(1-daf)+0.5f;
			
//			refWeights[i]=1;
//			insWeights[i]=1;
//			delWeights[i]=1;
		}
//		assert(false) : refWeights[0];
	}
	
	//TODO:  This does not seem to work properly; debug it.
	public int findBestOrientation(Read r, float minRatio){
		FlatAligner2 fla=new FlatAligner2();
		AlignmentResult ar0=fla.alignForwardShort(r.bases, original, 0, original.length-1, minRatio);
		AminoAcid.reverseComplementBasesInPlace(r.bases);
		AlignmentResult ar1=fla.alignForwardShort(r.bases, original, 0, original.length-1, minRatio);
		AminoAcid.reverseComplementBasesInPlace(r.bases);
		AlignmentResult ar=(ar0==null ? ar1 : ar1==null ? ar0 : ar0.ratio>=ar1.ratio ? ar0 : ar1);
		return ar==null ? -1 : ar==ar0 ? 0 : 1;
	}
	
	/** Generates match string and returns identity */
	public float alignAndGenerateMatch(Read r, Aligner ssa){
		return alignAndGenerateMatch(r, ssa, false, -1);
	}
	
	/** Generates match string and returns identity */
	public float alignAndGenerateMatch(Read r, Aligner ssa, boolean findOrientation, float minRatio){
		if(findOrientation){
			int orientation=findBestOrientation(r, minRatio);
			if(orientation<0){return 0;}
			else if(orientation==1){
				r.reverseComplement();
				r.setStrand(r.strand()^1);
			}
		}
		if(ssa==null){ssa=GeneCaller.getSSA();}
		byte[] query=r.bases;
		int a=0, b=ref.length-1;
		int[] max=ssa.fillUnlimited(query, original, a, b, -9999);
		if(max==null){return 0;}
		
		final int rows=max[0];
		final int maxCol=max[1];
		final int maxState=max[2];
		final byte[] match=ssa.traceback(query, original, a, b, rows, maxCol, maxState);
		int[] score=ssa.score(query, original, a, b, rows, maxCol, 0);
		r.match=match;
		r.start=score[1];
		r.stop=score[2];
		r.setMapped(true);
		final float identity=Read.identity(match);
		return identity;
	}
	
	void calcProbs(){
		for(BaseNode n : ref){
			for(BaseNode i=n; i!=null; i=i.insEdge){
				i.calcProbs();
			}
		}
	}
	
//	void makeWeights(){
////		assert(insWeights==null);
////		refWeights=new float[ref.length];
//		insWeights=new float[ref.length];
//		delWeights=new float[ref.length];
////		Arrays.fill(refWeights, 1f);
//		Arrays.fill(insWeights, 1f);
//		Arrays.fill(delWeights, 1f);
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public long readTotal;
	public long baseTotal;
	public long symbolTotal;
	public long msTotal;
	public long insTotal;
	public long delTotal;
	public long qualTotal;
	
	
	public int subCount=0;
	public int refCount=0;
	public int delCount=0;
	public int insCount=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public final String name;
	public final byte[] original;
	public final long numericID;
	public final int pad;
	//For ref nodes, calculate total outgoing weight.
	//Choose ins only if it is plurality allele. (or alternatively, majority allele).
	//But return to ref once it is no longer the majority/plurality of the outgoing weight.
	
	public final BaseNode[] ref;
	public final BaseNode[] del;

	//Original weights corresponding to conservation of sequence
	//This is for alignment, not scoring or updating, and is a simplified version of the graph
	//based on quality scores.
	public float[] refWeights;
//	public float[] insWeights;
//	public float[] delWeights;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final BaseNode dummy=new BaseNode('.', INS, 0);
	public static byte fakeQuality=15;
	public static boolean qualityByMS=true;
	
}
