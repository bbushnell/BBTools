package bin;

import aligner.IDAligner;
import fileIO.ByteStreamWriter;
import ml.CellNet;
import shared.Tools;
import shared.Vector;
import structures.ByteBuilder;
import structures.FloatList;
import tax.TaxTree;

public class Oracle extends BinObject implements Cloneable {

	public Oracle(float maxGCDif_, float maxDepthRatio_, float max3merDif_, float max4merDif_, float max5merDif_,
			float maxProduct_, float maxCovariance_, float minKmerProb_, int minEdgeWeight_) {
		maxGCDif0=maxGCDif_;
		maxDepthRatio0=maxDepthRatio_;
		max3merDif0=max3merDif_;
		max4merDif0=max4merDif_;
		max5merDif0=max5merDif_;
		maxProduct0=maxProduct_;
		maxCovariance0=maxCovariance_;
		minEdgeWeight=minEdgeWeight_;
		stringency0=1;
		minKmerProb0=minKmerProb_;
		if(BinObject.net0small!=null) {networkSmall=BinObject.net0small.copy(false);}
		if(BinObject.net0mid!=null) {networkMid=BinObject.net0mid.copy(false);}
		if(BinObject.net0large!=null) {networkLarge=BinObject.net0large.copy(false);}
	}
	
	public Oracle(float stringency, int minEdgeWeight_) {
		max3merDif0=Binner.max3merDif2*stringency;
		max4merDif0=Binner.max4merDif2*stringency;
		max5merDif0=Binner.max5merDif2*stringency;
		maxDepthRatio0=1+((Binner.maxDepthRatio2-1)*stringency);
		maxGCDif0=Binner.maxGCDif2*stringency;
		maxProduct0=max4merDif0*maxDepthRatio0*Binner.productMult;
		maxCovariance0=Binner.maxCovariance2*stringency;
		minKmerProb0=Binner.minKmerProb2;
		minEdgeWeight=minEdgeWeight_;
		stringency0=stringency;
		if(BinObject.net0small!=null) {networkSmall=BinObject.net0small.copy(false);}
		if(BinObject.net0mid!=null) {networkMid=BinObject.net0mid.copy(false);}
		if(BinObject.net0large!=null) {networkLarge=BinObject.net0large.copy(false);}
	}
	
	void clear() {
		best=null;
		bestIdx=-1;
		score=topScore=-1;
	}
	
	/** Higher is more similar */
	static final float similarity(float ratio_, float gcDif_, 
			float simDif_, float covariance_, float kmerProb_, long edges_) {
		final float ratio=ratio_;
		final float gcDif=gcDif_+1f;
		final float simDif=simDif_*0.5f+1f;
		final float covariance=1+covariance_*32;
		float product=simDif*ratio*gcDif*covariance/kmerProb_;
		if(BinObject.verbose) {
			System.err.println(product+"="+simDif+"*"+ratio+"*"+gcDif+"*"+covariance+"/"+kmerProb_);
		}
		return 1f/product;
	}
	
	public final float similarity(Bin a, Bin b, float stringency0) {
		long size=Tools.min(a.size(), b.size());
		float sizeMult=Binner.sizeAdjustMult(size);
		float stringency=stringency0*sizeMult;
		if(a.hasSSU() && b.hasSSU()) {stringency*=2;}

		if(a.maxDepth()==0 || b.maxDepth()==0) {
			stringency*=0.8f;//Has no effect...?  Maybe it will on coassemblies though.
		}
		float max3merDif=max3merDif0*stringency;
		float max4merDif=max4merDif0*stringency;
		float max5merDif=max5merDif0*stringency;
		float maxDepthRatio=1+((maxDepthRatio0-1)*stringency);
		float maxGCDif=maxGCDif0*stringency;
		float maxProduct=max4merDif*maxDepthRatio*Binner.productMult;
		maxProduct=Math.max(maxProduct, maxProduct0*sizeMult);
		float maxCovariance=maxCovariance0*stringency;
		float minProb=1-stringency0*(1-minKmerProb0);
		score=similarity(a, b, maxGCDif, maxDepthRatio, max3merDif, max4merDif, max5merDif,
				maxProduct, maxCovariance, minProb);
		return score;
	}
	
	public static float edgeMult(long e1, long e2, long eT, float d1, float d2) {
		long minEdges=Tools.min(e1, e2);
		if(minEdges<Binner.minEdgeWeight) {return eT<Binner.minEdgeWeight ? 1f : Binner.goodTransEdgeMult;}
		long maxEdges=Tools.max(e1, e2);
		float depth=0.5f*(d1+d2);
		if(minEdges>Binner.lowDepthEdgeRatio*depth && maxEdges<Binner.highDepthEdgeRatio*depth) {
			return Binner.goodEdgeMult;
		}
		return eT<Binner.minEdgeWeight ? 1f : Binner.goodTransEdgeMult;
	}
	
	/** Higher is more similar */
	private final float similarity(Bin a, Bin b, float maxGCDif, float maxDepthRatio, float max3merDif, 
			float max4merDif, float max5merDif, float maxProduct, float maxCovariance, float minKmerProb) {
//		if(bsw!=null) {return similarityB(a, b, maxGCDif, maxDepthRatio, maxKmerDif, maxProduct, maxCovariance);}
		fastComparisons++;
		
		if(BinObject.verbose || verbose2) {
			System.err.println("Comparing to "+b.id()+": "+
					"maxKmerDif="+max4merDif+", maxDepthRatio="+maxDepthRatio+
					", maxProduct="+maxProduct+", maxGCDif="+maxGCDif+
					", maxCovariance="+maxCovariance);
		}
		final boolean sameLabel=(a.labelTaxid>0 && a.labelTaxid==b.labelTaxid);
//		final boolean diffLabel=(a.labelTaxid<1 || a.labelTaxid!=b.labelTaxid);
		if(Binner.PERFECT_ORACLE) {return sameLabel ? 1-1f/b.size() : -1;}
		
		final float gcDif=Math.abs(a.gc()-b.gc());
		final float depthRatio=a.depthRatio(b);
		final long minlen=Math.min(a.size(), b.size());
		if(BinObject.verbose || verbose2) {
			System.err.println("gcdif="+gcDif);
			System.err.println("depthRatio="+depthRatio);
		}
		if(gcDif>maxGCDif*Binner.goodEdgeMult || 
				depthRatio>maxDepthRatio*Binner.goodEdgeMult) {
			return -1;
		}//Early exit before edge-tracking
		
		final long edges1=useEdges ? a.countEdgesTo(b) : 0;
		final long edges2=useEdges ? b.countEdgesTo(a) : 0;
		final long minEdges=Tools.min(edges1, edges2);
		if(BinObject.verbose || verbose2) {
			System.err.println("A: size="+minlen+", e1="+edges1+", e2="+edges2+", minWeight="+minEdgeWeight);
		}
		final long edgesT=(minEdges>=Binner.minEdgeWeight ? minEdges : a.transEdgesTo(b));
//		if(minEdges<1 && minlen<3000) {
//			return -1;
//		}//Doesn't help or even do much at minlen=3000
		float mult=edgeMult(edges1, edges2, edgesT, a.depthTotal(), b.depthTotal());
//		mult*=Math.sqrt(Binner.sizeAdjustMult(minlen));
//		if(minlen<3000) {mult*=0.0f;return -1;}
		
		if(BinObject.verbose || verbose2) {
			System.err.println("B: mult="+mult+", gcdif="+gcDif+", max="+(maxGCDif*mult));
		}
		if(gcDif>maxGCDif*mult*Binner.cutoffMultD) {
//			assert(!sameLabel) : "gcdif="+gcDif+">"+(maxGCDif*mult*Binner.cutoffMultD);
			return -1;
		}
		float covariance=a.covariance(b);
		if(BinObject.verbose || verbose2) {
			System.err.println("C: depthRatio="+depthRatio+", max="+(maxDepthRatio*mult)+
					", covariance="+covariance+", max="+(maxCovariance*mult));
		}
		if(depthRatio>maxDepthRatio*mult*Binner.cutoffMultD || 
				covariance>maxCovariance*mult*Binner.cutoffMultD) {return -1;}
//		if(!taxaOK(a.taxid(), b.taxid())) {return -1;}

		trimerComparisons++;
		float trimerDif=(countTrimers ? 
				Vector.cosineDifference(a.trimers, b.trimers) : 0);
//				SimilarityMeasures.cosineDifference(a.trimers, b.trimers) : 0);
		
		//This causes a large speedup by avoiding tetramer calculation
		//0.75 has no effect, so 0.8 is safe (0.725 causes a slight change) 
		if(trimerDif>max3merDif*mult*Binner.cutoffMultA ||
				trimerDif*depthRatio>maxProduct*mult*Binner.cutoffMultB*0.8f) {return -1;}

		tetramerComparisons++;
		float tetramerDif=Vector.cosineDifference(a.tetramers, b.tetramers);
		final float product=tetramerDif*depthRatio;
		float kmerProb=KmerProb.prob(minlen, tetramerDif);
		kmerProb=1-(1-kmerProb)/mult;
		if(tetramerDif>max4merDif*mult*Binner.cutoffMultA || 
				product>maxProduct*mult*Binner.cutoffMultB || kmerProb<0.5f) {return -1;}

		slowComparisons++;
		float pentamerDif=(a.numPentamers<BinObject.minPentamerSizeCompare ||
				b.numPentamers<BinObject.minPentamerSizeCompare ? Math.min(1, tetramerDif*1.7f) :
					Vector.cosineDifference(a.pentamers, b.pentamers));
		if(BinObject.verbose || verbose2) {
			System.err.println("D: tetramerDif="+tetramerDif+", max="+(max4merDif*mult)+
					", product="+product+", max="+(maxProduct*mult));
		}
		if(pentamerDif>max5merDif*mult*Binner.cutoffMultA) {return -1;}
		
		final float similarity=similarity(depthRatio, gcDif, tetramerDif, covariance, kmerProb, 
				Tools.min(edges1, edges2));
		final float netOutput=runNetwork(a, b, minEdges, edgesT, gcDif, depthRatio, covariance, 
				trimerDif, tetramerDif, pentamerDif, kmerProb, similarity, false);
		
		if(BinObject.verbose || verbose2) {
			System.err.println("E: similarity="+similarity+", netOutput="+netOutput);
		}
		
		CellNet network=getNetwork(minlen);
		float mult2=mult;
		if(network!=null) {
			final float cutoff=(makingBinMap && Binner.netCutoff1>network.cutoff ? 
					Binner.netCutoff1 : network.cutoff);
			netComparisons++;
			
			if(netOutput>Binner.netCutoffUpper) {
				mult2*=Binner.netMultUpper;
			}else if(netOutput<Binner.netCutoffLower) {
				mult2*=Binner.netMultLower;
			}
			float ratio=(netOutput<0.001f ? 0 : netOutput/cutoff);
//			mult2=(float)(mult2*ratio*ratio*Math.sqrt(ratio));
			mult2=(float)(mult2*ratio*ratio);
		}
		if(BinObject.verbose || verbose2) {
			System.err.println("F: mult="+mult+", tetramerDif="+tetramerDif+", max="+(max4merDif*mult)+
					", product="+product+", max="+(maxProduct*mult)+
					"\n   kmerProb="+kmerProb+", min="+(minKmerProb));
		}
		
		float ret=netOutput;
		if(trimerDif>max3merDif*mult2 || tetramerDif>max4merDif*mult2 || pentamerDif>max5merDif*mult2 || 
				product>maxProduct*mult2 || kmerProb<minKmerProb) {ret=-1;}

		float mult3=(network==null ? mult : mult2*Binner.cutoffMultC);
		if(gcDif>maxGCDif*mult3 || depthRatio>maxDepthRatio*mult3 || covariance>maxCovariance*mult3) {
			ret=-1;
		}
		
		if(bsw!=null && canEmitVector(a, b, ret)) {
			if(sameLabel || Math.random()<=negativeEmitProb) {emitVector(a, b, bsw);}
		}
		
		if(Binner.BAN_BAD_MERGES && !sameLabel) {ret=-1;}
		if(ret>-1 && ssa!=null) {
			float id=ssuCompatibility(a, b);
			if(id<minSSUID) {ret=-1;}
			else if(id<2) {ret+=id;}
		}
		return ret;
	}
	
	final float ssuCompatibility(Bin a, Bin b) {
		if(a.r16S==null && a.r18S==null) {return 2;}
		if(b.r16S==null && b.r18S==null) {return 2;}
		if(a.r16S!=null && b.r18S!=null) {return -1;}
		if(a.r18S!=null && b.r16S!=null) {return -1;}
		if(a.r16S!=null && b.r16S!=null) {return ssa.align(a.r16S, b.r16S, null, 0);}
		if(a.r18S!=null && b.r18S!=null) {return ssa.align(a.r18S, b.r18S, null, 0);}
		assert(false);
		return 0;
	}
	
	final float runNetwork(Bin a, Bin b, final long minEdges, final long transEdges, final float gcDif, 
			final float depthRatio, final float covariance, final float trimerDif, final float tetramerDif,
			final float pentamerDif, final float kmerProb, final float similarity, boolean includeAnswer) {
		CellNet network=getNetwork(Math.min(a.size(), b.size()));
		if(network==null) {return similarity;}
		if(vector==null) {vector=new FloatList();}
		toVector(a, b, minEdges, transEdges, gcDif, depthRatio,
				covariance, trimerDif, tetramerDif, pentamerDif, kmerProb, similarity, vector, false);
		network.applyInput(vector);
		float result=network.feedForward();
		return result;
	}
	
	static String header() {
		ByteBuilder bb=new ByteBuilder();
//		bb.append('#').append("aSize").tab().append("bSize");//0 1
//		bb.tab().append("aGC").tab().append("bGC");//2 3
//		bb.tab().append("gcDif");//4
//		bb.tab().append("depthRatio").tab().append("covariance");//5 6
//		bb.tab().append("aDepth").tab().append("bDepth");//7 8
//		bb.tab().append("numDepth").tab().append("tetramerDif");//9 10
//		bb.tab().append("aEntrop").tab().append("bEntrop");//11 12
//		bb.tab().append("entDif");//13
//		bb.tab().append("aSpec").tab().append("bSpec");//14 15
//		bb.tab().append("specDif");//16
//		bb.tab().append("aContigs").tab().append("bContigs");//17 18
//		bb.tab().append("aEdgeW").tab().append("bEdgeW");//19 20
//		bb.tab().append("aEdges").tab().append("bEdges");//21 22
//		bb.tab().append("similarity");//23
//		bb.tab().append("sameTax");//24
		

		bb.append('#').append("aSize");//0
		bb.tab().append("bGC");//1
		bb.tab().append("gcDif");//2
		bb.tab().append("depthRatio").tab().append("covariance");//3 4
		bb.tab().append("bDepth");//5
		bb.tab().append("numDepth").tab().append("tetramerDif");//6 7
		bb.tab().append("bEntrop");//8
		bb.tab().append("entDif");//9
		bb.tab().append("bSpec");//10
		bb.tab().append("specDif");//11
		bb.tab().append("minEdge");//12
		bb.tab().append("similarity");//13
		bb.tab().append("sameTax");//14
		
		return bb.toString();
	}
	
	FloatList toVector(Bin a, Bin b, FloatList list, boolean includeAnswer) {
		if(a.size()>b.size()) {return toVector(b, a, list, includeAnswer);}
		if(list==null) {list=new FloatList();}
		list.clear();
		long minlen=Math.min(a.size(), b.size());
		long edges1=a.countEdgesTo(b);
		long edges2=b.countEdgesTo(a);
		long minEdges=Tools.min(edges1, edges2);
		final long edgesT=(minEdges>=Binner.minEdgeWeight ? minEdges : a.transEdgesTo(b));
		float gcDif=Math.abs(a.gc()-b.gc());
		float depthRatio=a.depthRatio(b);
		float covariance=a.covariance(b);
		float tetramerDif=SimilarityMeasures.cosineDifference(a.tetramers, b.tetramers);
		float trimerDif=(countTrimers ? 
				SimilarityMeasures.cosineDifference(a.trimers, b.trimers) : 0f);
		float pentamerDif=(a.numPentamers<BinObject.minPentamerSizeCompare ||
				b.numPentamers<BinObject.minPentamerSizeCompare ? Math.min(1, tetramerDif*1.7f) :
					SimilarityMeasures.cosineDifference(a.pentamers, b.pentamers));
		float kmerProb=KmerProb.prob(minlen, tetramerDif);
		float similarity=similarity(depthRatio, gcDif, tetramerDif, covariance, kmerProb, Tools.min(edges1, edges2));
		
		return toVector(a, b, minEdges, edgesT, gcDif, depthRatio, covariance, 
				trimerDif, tetramerDif, pentamerDif, kmerProb, similarity, list, includeAnswer);
	}
	
	FloatList toVector(Bin a, Bin b, final long minEdges, final long transEdges, final float gcDif, final float depthRatio,
			final float covariance, float trimerDif, float tetramerDif, float pentamerDif, float kmerProb, 
			final float similarity, FloatList list, boolean includeAnswer) {
		if(a.size()>b.size()) {
			return toVector(b, a, minEdges, transEdges, gcDif, depthRatio, covariance, 
					trimerDif, tetramerDif, pentamerDif, kmerProb, similarity, list, includeAnswer);
		}
		if(list==null) {list=new FloatList();}
		list.clear();
		float depth=0.5f*(a.depthTotal()+b.depthTotal()+0.5f);
		
		int numDepths=BinObject.samplesEquivalent;
//		long minlen=Math.min(a.size(), b.size());
//		final float kmerProb=KmerProb.prob(minlen, tetramerDif);
		float logSize=(float)Tools.log2(Tools.max(128, a.size()));
		float logSize2=(float)Tools.log2(Tools.max(128, b.size()));
		float sizeProxy=0.1f*(logSize-7);
		float sizeProxy2=0.1f*(logSize2-7);
		
		float depthRatioProxy=(float)(0.5f*Tools.log2(depthRatio));
		if(depthRatioMethod!=1) {depthRatioProxy=a.depthRatio(b, depthRatioMethod);}
//		if(depthRatioProxy<=0) {
//			KillSwitch.kill(depthRatioMethod+", "+
//					depthRatio+", "+depthRatioProxy+", "+a.depth(0)+", "+b.depth(0));
//		}
		final boolean pentamers=(a.numPentamers>=BinObject.minPentamerSizeCompare &&
				b.numPentamers<BinObject.minPentamerSizeCompare);
		
//		float euc=(addEuclidian ? SimilarityMeasures.euclideanDistance(a.tetramers, b.tetramers) : 0);
		final float invA4=Tools.invSum(a.tetramers), invB4=Tools.invSum(b.tetramers);
		float hel=(addHellinger ? SimilarityMeasures.hellingerDistance(a.tetramers, b.tetramers, invA4, invB4) : 0);
		float hel3=(addHellinger3 ? SimilarityMeasures.hellingerDistance(a.trimers, b.trimers) : 0);
		float hel5=(!addHellinger5 ? 0 : !pentamers ? Math.min(1, hel*1.7f) :
					SimilarityMeasures.hellingerDistance(a.pentamers, b.pentamers));
		float jsdiv=(addJsDiv ? SimilarityMeasures.jensenShannonDivergence(a.tetramers, b.tetramers, invA4, invB4) : 0);
		float absdif=(addAbsDif ? SimilarityMeasures.absDif(a.tetramers, b.tetramers, invA4, invB4) : 0);
		float mult=vectorSmallNumberMult;//For making very small numbers bigger
		float gcComp=(!addGCComp ? 0 : 
			SimilarityMeasures.cosineDifferenceCompensated(a.tetramers, b.tetramers, 4));
		
		if(vectorSmallNumberRoot) {
			trimerDif=(float)Math.sqrt(trimerDif);
			tetramerDif=(float)Math.sqrt(tetramerDif);
			pentamerDif=(float)Math.sqrt(pentamerDif);
			absdif=(float)Math.sqrt(absdif);
			jsdiv=(float)Math.sqrt(jsdiv);
			hel=(float)Math.sqrt(hel);
			hel3=(float)Math.sqrt(hel3);
			hel5=(float)Math.sqrt(hel5);
		}
		
//		list.add(a.size());//-1
		if(printSizeInVector) {
			list.add(a.size());
			list.add(a.numContigs());
			list.add(b.size());
			list.add(b.numContigs());
		}
		list.add(sizeProxy);//0
		list.add(sizeProxy2);//1
//		list.add(0.1f*(float)Tools.log2(b.size()));
//		list.add(a.gc());
		list.add(b.gc());//2
		list.add(gcDif);//3

		list.add(trimerDif*mult);//4
		list.add(tetramerDif*mult);//5
		list.add(pentamerDif*mult);//6
		list.add(absdif*mult);//7
		list.add(jsdiv*mult);//8
		list.add(hel*mult);//9
		list.add(hel3*mult);//10
		list.add(hel5*mult);//11
		list.add(gcComp*mult);//12
		
		list.add(kmerProb);//13
		list.add((float)(minEdges/depth));//14
		list.add((float)(transEdges/depth));//15
		list.add(depthRatioProxy);//16
		list.add((float)Math.sqrt(covariance)*mult);//17
//		list.add((float)(0.1f*Tools.log2(1+a.depthTotal())));
		list.add((float)(0.1f*Tools.log2(1+b.depthTotal())));//18
		
		list.add(numDepths>1 ? 1 : 0);//19
		list.add(numDepths>2 ? 1 : 0);//20
		list.add(numDepths>3 ? 1 : 0);//21
		list.add((float)(1.25f*(sampleEntropy-1)/(sampleEntropy+2)));//22 //sampleEntropy is like numDepths
//		list.add(0.2f*((float)sampleEntropy-1));//15
		
//		list.add(a.entropy);
		list.add(addEntropy ? 8*(1-b.entropy) : 0);//23
		list.add(addEntropy ? 8*Tools.absdif(a.entropy, b.entropy) : 0);//24
//		list.add(a.strandedness-1);
		list.add(addStrandedness ? b.strandedness : 0);//25
		list.add(addStrandedness ? 8*Tools.absdif(a.strandedness, b.strandedness) : 0);//26
//		list.add((float)(0.5f*Tools.log2(a.numContigs())));
//		list.add((float)(0.5f*Tools.log2(b.numContigs())));
//		list.add(edges1/depth);
//		list.add(edges2/depth);
//		list.add(0.1f*a.numEdges()/(float)(Tools.max(b.numContigs(), 1)));
//		list.add(0.1f*b.numEdges()/(float)(Tools.max(b.numContigs(), 1)));
		list.add(1-similarity);//27
		if(printNetOutputInVector) {
			float out1=0, out2=0, out3=0;
			if(networkSmall!=null) {
				networkSmall.applyInput(list);
				out1=networkSmall.feedForward();
			}
			if(networkMid!=null && networkMid!=networkSmall) {
				networkMid.applyInput(list);
				out2=networkMid.feedForward();
			}
			if(networkLarge!=null && networkLarge!=networkMid) {
				networkLarge.applyInput(list);
				out3=networkLarge.feedForward();
			}
			list.add(out1);
			list.add(out2);
			list.add(out3);
		}
		if(includeAnswer) {
			if(printWeightInVector>0) {
				if(printWeightInVector==1) {
					list.add(0.125f+sizeProxy);
				}else if(printWeightInVector==2) {
					list.add(0.125f+sizeProxy*sizeProxy);
				}else if(printWeightInVector==3) {
					list.add(0.125f+(float)Math.sqrt(a.size()/20000f));
				}else {
					assert(false) : printWeightInVector;
				}
			}
			assert(a.labelTaxid>0 && b.labelTaxid>0) : a.labelTaxid+", "+b.labelTaxid+", "+a.name();
			list.add(a.labelTaxid==b.labelTaxid ? 1 : 0);
		}
		list.shrink();
		for(float f : list.array) {
			assert(Float.isFinite(f)) : list;
		}
		return list;
	}
	
	private boolean taxaOK(int aTaxid, int bTaxid) {
		if(!allowHalfTaxID && (aTaxid<1 || bTaxid<1)) {return false;}
		if(!allowNoTaxID && aTaxid<1 && bTaxid<1) {return false;}
		if(taxlevel<0 || BinObject.tree==null || aTaxid==bTaxid || aTaxid<1 || bTaxid<1) {return true;}
		int commonAncestorLevel=BinObject.tree.commonAncestorLevel(aTaxid, bTaxid);
		return (commonAncestorLevel<=taxlevel);
	}
	
	
	protected Oracle clone() {
		try {
			Oracle clone=(Oracle) super.clone();
			clone.best=null;
			clone.vector=null;
			clone.networkSmall=(networkSmall==null ? null : networkSmall.copy(false));
			clone.networkMid=(networkMid==null ? null : networkMid.copy(false));
			clone.networkLarge=(networkLarge==null ? null : networkLarge.copy(false));
			clone.ssa=(ssa==null ? null : aligner.Factory.makeIDAligner());
			return clone;
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}
	
	static boolean canEmitVector(Bin a, Bin b, float result) {
		long minSize=Math.min(a.size(), b.size());
		if(a.labelTaxid<1 || b.labelTaxid<1) {return false;}
		if(minSize<minEmitSize || minSize>maxEmitSize) {return false;}
		if(!a.pure() || !b.pure()) {return false;}
		boolean sameTID=(a.labelTaxid==b.labelTaxid);
		boolean merge=(result>0);
		final boolean ret;
		if(merge) {ret=(sameTID ? emitTP : emitFP);}
		else{ret=(sameTID ? (emitFN && !makingBinMap) : emitFP);}
//		assert(!ret) : result+", "+merge+", "+sameTID+", "+minSize;
		return ret;
	}
	
	void emitVector(Bin a, Bin b, ByteStreamWriter bsw) {
		assert(bsw!=null);
		long minlen=Tools.min(a.size(), b.size());
		float prob=(minlen/4000f)*(minlen/2000f);
		if(Math.random()>prob) {return;}
		if(vector==null) {vector=new FloatList();}
		toVector(a, b, vector, true);
		emitVector(vector, bsw);
	}
	
	static void emitVector(FloatList vector, ByteStreamWriter bsw) {
		assert(bsw!=null);
		ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<vector.size(); i++) {
			if(i>0) {bb.tab();}
			bb.append(vector.get(i), 7, true);
		}
		synchronized(bsw) {
			bsw.print(bb.nl());
		}
	}
	
	private CellNet getNetwork(long size) {
		if(size<Binner.minNetSize) {return null;}
		if(size<Binner.midNetSize) {return networkSmall;}
		if(size<Binner.largeNetSize) {return networkMid;}
		return networkLarge;
	}
	
	Bin best=null;
	float score=-1;
	float topScore=-1;
	int bestIdx=-1;
	
	long fastComparisons=0;
	long trimerComparisons=0;
	long tetramerComparisons=0;
	long slowComparisons=0;
	long netComparisons=0;

	final float max3merDif0;
	final float max4merDif0;
	final float max5merDif0;
	final float maxDepthRatio0;
	final float maxGCDif0;
	final float maxProduct0;
	final float maxCovariance0;
	final float minKmerProb0;
	final int minEdgeWeight;
	
	final float stringency0;
	private FloatList vector;
	private CellNet networkSmall;
	private CellNet networkMid;
	private CellNet networkLarge;
	private IDAligner ssa=(SpectraCounter.call16S ? 
			aligner.Factory.makeIDAligner() : null);
	
	int taxlevel=TaxTree.SPECIES;
	boolean allowNoTaxID=true;
	boolean allowHalfTaxID=true;
	boolean useEdges=true;
	
	static ByteStreamWriter bsw;
	static boolean emitTP=true;
	static boolean emitFP=true;
	static boolean emitTN=true;
	static boolean emitFN=true;
	static int minEmitSize=0;
	static int maxEmitSize=2000000000;
	static double negativeEmitProb=1;
	static boolean printSizeInVector=false;
	static int printWeightInVector=1;
	static boolean printNetOutputInVector=false;
	static float minSSUID=0.98f;
	boolean verbose2=false;
	
	
}
