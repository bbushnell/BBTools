package align2;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Aug 5, 2013
 *
 */
public class BandedAlignerConcrete extends BandedAligner{

	//The live banded edit-distance aligner (makeBandedAligner always returns this). Two rolling rows
	//(array1/array2 swapped via arrayTemp) give O(width) memory; each cell = min(up, diag, left) with
	//+1 per mismatch/indel and 0 for a match (forceDiag on the last row and at the ref boundary).
	//The four methods are mirror images — forward/reverse flip column iteration and the force-diagonal
	//boundary (col==ref.length-1 vs col==0); the RC variants complement the query base. All debug prints
	//are verbose-guarded. Bounds + mirror-symmetry independently adversarially verified clean (2026-06-20).

	/**
	 * Testing harness for banded alignment algorithms.
	 * Runs forward, reverse, and combination alignment tests with both penalized
	 * and unpenalized off-center scoring modes.
	 * @param args Command line arguments: query_seq ref_seq [qstart] [rstart] [maxedits] [width]
	 */
	public static void main(String[] args){
		byte[] query=args[0].getBytes();
		byte[] ref=(args[1].equals(".") ? args[0].getBytes() : args[1].getBytes());
		int qstart=-1;
		int rstart=-1;
		int maxedits=big-1;
		int width=5;
		if(args.length>2){qstart=Integer.parseInt(args[2]);}
		if(args.length>3){rstart=Integer.parseInt(args[3]);}
		if(args.length>4){maxedits=Integer.parseInt(args[4]);}
		if(args.length>5){width=Integer.parseInt(args[5]);}
		
		BandedAlignerConcrete ba=new BandedAlignerConcrete(width);
		
		int edits;
		
		penalizeOffCenter=true;
		edits=ba.alignForward(query, ref, (qstart==-1 ? 0 : qstart), (rstart==-1 ? 0 : rstart), maxedits, true);
		System.out.println("Forward:    \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
		System.out.println("***********************\n");

		penalizeOffCenter=false;
		edits=ba.alignForward(query, ref, (qstart==-1 ? 0 : qstart), (rstart==-1 ? 0 : rstart), maxedits, true);
		System.out.println("Forward2:   \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
		System.out.println("***********************\n");
//
//		edits=ba.alignForwardRC(query, ref, (qstart==-1 ? query.length-1 : qstart), (rstart==-1 ? 0 : rstart), maxedits, true);
//		System.out.println("ForwardRC:  \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
//		System.out.println("***********************\n");

		penalizeOffCenter=true;
		edits=ba.alignReverse(query, ref, (qstart==-1 ? query.length-1 : qstart), (rstart==-1 ? ref.length-1 : rstart), maxedits, true);
		System.out.println("Reverse:    \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
		System.out.println("***********************\n");

		penalizeOffCenter=false;
		edits=ba.alignReverse(query, ref, (qstart==-1 ? query.length-1 : qstart), (rstart==-1 ? ref.length-1 : rstart), maxedits, true);
		System.out.println("Reverse2:   \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
		System.out.println("***********************\n");
		
//		edits=ba.alignReverseRC(query, ref, (qstart==-1 ? 0 : qstart), (rstart==-1 ? ref.length-1 : rstart), maxedits, true);
//		System.out.println("ReverseRC:  \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
//		System.out.println("***********************\n");

		penalizeOffCenter=true;
		edits=ba.alignQuadruple(query, ref, maxedits, true);
		System.out.println("Quadruple:    \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
		System.out.println("***********************\n");

		penalizeOffCenter=false;
		edits=ba.alignQuadruple(query, ref, maxedits, true);
		System.out.println("Quadruple2:   \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
		System.out.println("***********************\n");

		penalizeOffCenter=true;
		edits=ba.alignDouble(query, ref, maxedits, true);
		System.out.println("Double:    \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
		System.out.println("***********************\n");

		penalizeOffCenter=false;
		edits=ba.alignDouble(query, ref, maxedits, true);
		System.out.println("Double2:    \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
		System.out.println("***********************\n");
	}
	
	
	/**
	 * Creates a new banded aligner with specified band width.
	 * Initializes working arrays and validates that the band width is sufficient
	 * for the maximum edit distance calculations.
	 * @param width_ Band width for alignment (must be odd, minimum 3)
	 */
	public BandedAlignerConcrete(int width_){
		super(width_);
		array1=new int[maxWidth+2];
		array2=new int[maxWidth+2];
		//Size maxWidth+2 IS the bounds guarantee: per-row width<=maxWidth, and the DP indexes mloc in
		//[1, width], touching mloc-1>=0 and mloc+1<=width+1<=maxWidth+1 (last valid index). The whole-
		//array big fill here is load-bearing: the band-edge read arrayPrev[width+1] is NEVER re-written
		//by the per-row loop, so it always returns this 'big' sentinel ("no predecessor"). Verified.
		Arrays.fill(array1, big);
		Arrays.fill(array2, big);
//		for(int i=2; i<rows; i++){
//			matrix[i]=matrix[i-2];
//		}
		assert(big>maxWidth/2);
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	@Override
	public int alignForward(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact){
		assert(big>maxEdits);
		if(verbose){System.err.println("alignForward("+new String(query)+", "+new String(ref)+", "+qstart+", "+rstart+", "+maxEdits+")");}
		if(query.length-qstart>ref.length-rstart){
			int x=alignForward(ref, query, rstart, qstart, maxEdits, exact);
			int temp=lastQueryLoc;
			lastQueryLoc=lastRefLoc;
			lastRefLoc=temp;
			if(verbose){
				System.out.println("Reversed.");
				System.out.println("Final state: lastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
						", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+"\n");
			}
			return x;
		}
		int edits=0, row=0;
		lastRow=-1;
		lastEdits=0;
		lastOffset=0;
		
		final int width=Tools.min(maxWidth, (maxEdits*2)+1, Tools.max(query.length, ref.length)*2+2)|1;
		final int halfWidth=width/2;
		final boolean inexact=!exact;
		
		int qloc=qstart;
		int rsloc=rstart-halfWidth;
		final int xlines=query.length-qstart;
		final int ylines=ref.length-rstart;
		final int len=Tools.min(xlines, ylines);
		if(verbose){System.err.println("xlines="+xlines+", ylines="+ylines+", len="+len);}
		if(len<1){
			if(false){
				throw new RuntimeException("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			}
			assert(false) : ("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			return 0;
		}

		Arrays.fill(array1, 0, Tools.min(width, maxWidth)+1, big);
		Arrays.fill(array2, 0, Tools.min(width, maxWidth)+1, big);
		arrayCurrent=array1;
		arrayPrev=array2;
		{
			if(verbose){System.err.println("\nFirst row.");}
			final byte q=query[qloc];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			edits=big;
			int mloc=1+(colStart-rsloc);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
//			assert(false) : mloc+", "+colStart+", "+rsloc;
			for(int col=colStart; col<colLimit; mloc++, col++){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int score=(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; prev=0; score="+score+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			row++; qloc++; rsloc++;
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}
		
		for(row=1; row<len; row++, qloc++, rsloc++){
//			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayCurrent));}
			arrayTemp=arrayCurrent;
			arrayCurrent=arrayPrev;
			arrayPrev=arrayTemp;
			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayPrev)+", current="+Arrays.toString(arrayCurrent));}
			final byte q=query[qloc];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			Arrays.fill(arrayCurrent, big);
			edits=big;
			int mloc=1+(colStart-rsloc);
			boolean forceDiag=(row==len-1);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colStart; col<colLimit; mloc++, col++){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int scoreUp=arrayPrev[mloc+1]+1;
				final int scoreDiag=arrayPrev[mloc]+(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				final int scoreLeft=arrayCurrent[mloc-1]+1;
				final int score=(forceDiag || col==ref.length-1) ? scoreDiag : Tools.min(scoreUp, scoreDiag, scoreLeft);
				if(verbose){System.err.println("prev=min(s["+(mloc-1)+"]="+arrayCurrent[mloc-1]+", p["+(mloc)+"]="+arrayPrev[mloc]+", p["+(mloc+1)+"]="+arrayPrev[mloc+1]+")");}
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; up="+scoreUp+"; diag="+scoreDiag+"; left="+scoreLeft+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			if(edits>maxEdits){row++; break;}
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}

		lastRow=row-1;
		lastEdits=edits;
		lastQueryLoc=qloc-1;
		lastOffset=lastOffset(arrayCurrent, halfWidth);
		lastRefLoc=rsloc+halfWidth-lastOffset-1;
		//Directional clamp (mirrored across the 4 methods): each guards only the overshoot direction
		//matching its traversal — forward clamps the high end, the reverse variants clamp the low end.
		//align2/BandedAlignerConcrete#001 (QUESTION, not a confirmed bug): the opposite direction is
		//structurally unguarded (could in principle exit with lastRefLoc<0), but the swap-guard geometry
		//(qlen-qstart<=rlen-rstart) is assumed to keep it unreachable. Adversarially verified: no bounds
		//violation here; lastRefLoc/lastQueryLoc are int fields, not array indices in this method.
		while(lastRefLoc>=ref.length || lastQueryLoc>=query.length){lastRefLoc--; lastQueryLoc--;}
		if(verbose){
			System.out.println("\nFinal state: arrayCurrent="+Arrays.toString(arrayCurrent)+"\nlastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
					", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+"\n");
		}
		return edits;
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	@Override
	public int alignForwardRC(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact){
		assert(big>maxEdits);
		if(verbose){System.err.println("alignForwardRC("+new String(query)+", "+new String(ref)+", "+qstart+", "+rstart+", "+maxEdits+")");}
		if(qstart+1>ref.length-rstart){
			int x=alignReverseRC(ref, query, rstart, qstart, maxEdits, exact);
			int temp=lastQueryLoc;
			lastQueryLoc=lastRefLoc;
			lastRefLoc=temp;
			if(verbose){
				System.out.println("Reversed.");
				System.out.println("Final state: lastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
						", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+"\n");
			}
			return x;
		}
		int edits=0, row=0;
		lastRow=-1;
		lastEdits=0;
		lastOffset=0;
		
		final int width=Tools.min(maxWidth, (maxEdits*2)+1, Tools.max(query.length, ref.length)*2+2)|1;
		final int halfWidth=width/2;
		final boolean inexact=!exact;
		
		int qloc=qstart;
		int rsloc=rstart-halfWidth;
		final int xlines=qstart+1;
		final int ylines=ref.length-rstart;
		final int len=Tools.min(xlines, ylines);
		if(verbose){System.err.println("xlines="+xlines+", ylines="+ylines+", len="+len);}
		if(len<1){
			if(false){
				throw new RuntimeException("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			}
			assert(false) : ("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			return 0;
		}

		Arrays.fill(array1, 0, Tools.min(width, maxWidth)+1, big);
		Arrays.fill(array2, 0, Tools.min(width, maxWidth)+1, big);
		arrayCurrent=array1;
		arrayPrev=array2;
		
		{
			if(verbose){System.err.println("\nFirst row.");}
			final byte q=AminoAcid.baseToComplementExtended[query[qloc]];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			edits=big;
			int mloc=1+(colStart-rsloc);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colStart; col<colLimit; mloc++, col++){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int score=(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; prev=0; score="+score+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			row++; qloc--; rsloc++;
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}
		
		for(row=1; row<len; row++, qloc--, rsloc++){
//			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayCurrent));}
			arrayTemp=arrayCurrent;
			arrayCurrent=arrayPrev;
			arrayPrev=arrayTemp;
			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayPrev)+", current="+Arrays.toString(arrayCurrent));}
			final byte q=AminoAcid.baseToComplementExtended[query[qloc]];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			Arrays.fill(arrayCurrent, big);
			edits=big;
			int mloc=1+(colStart-rsloc);
			boolean forceDiag=(row==len-1);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colStart; col<colLimit; mloc++, col++){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int scoreUp=arrayPrev[mloc+1]+1;
				final int scoreDiag=arrayPrev[mloc]+(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				final int scoreLeft=arrayCurrent[mloc-1]+1;
				final int score=(forceDiag || col==ref.length-1) ? scoreDiag : Tools.min(scoreUp, scoreDiag, scoreLeft);
				if(verbose){System.err.println("prev=min(s["+(mloc-1)+"]="+arrayCurrent[mloc-1]+", p["+(mloc)+"]="+arrayPrev[mloc]+", p["+(mloc+1)+"]="+arrayPrev[mloc+1]+")");}
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; up="+scoreUp+"; diag="+scoreDiag+"; left="+scoreLeft+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			if(edits>maxEdits){row++; break;}
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}

		lastRow=row-1;
		lastEdits=edits;
		lastOffset=lastOffset(arrayCurrent, halfWidth);
		lastQueryLoc=qloc+1;
		lastRefLoc=rsloc+halfWidth-lastOffset-1;
		while(lastRefLoc>=ref.length || lastQueryLoc<0){lastRefLoc--; lastQueryLoc++;}
		if(verbose){
			System.out.println("\nFinal state: arrayCurrent="+Arrays.toString(arrayCurrent)+"\nlastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
					", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+", qloc="+qloc+"\n");
		}
		return edits;
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	@Override
	public int alignReverse(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact){
		assert(big>maxEdits);
		if(verbose){System.err.println("alignReverse("+new String(query)+", "+new String(ref)+", "+qstart+", "+rstart+", "+maxEdits+")");}
		if(qstart>rstart){
			int x=alignReverse(ref, query, rstart, qstart, maxEdits, exact);
			int temp=lastQueryLoc;
			lastQueryLoc=lastRefLoc;
			lastRefLoc=temp;
			if(verbose){
				System.out.println("Reversed.");
				System.out.println("Final state: lastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
						", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+"\n");
			}
			return x;
		}
//		if(true){return big;}
		int edits=0, row=0;
		lastRow=-1;
		lastEdits=0;
		lastOffset=0;
		
		final int width=Tools.min(maxWidth, (maxEdits*2)+1, Tools.max(query.length, ref.length)*2+2)|1;
		final int halfWidth=width/2;
		final boolean inexact=!exact;
		
		int qloc=qstart;
		int rsloc=rstart-halfWidth;
		final int xlines=qstart+1;
		final int ylines=rstart+1;
		final int len=Tools.min(xlines, ylines);
		if(verbose){System.err.println("xlines="+xlines+", ylines="+ylines+", len="+len);}
		if(len<1){
			if(false){
				throw new RuntimeException("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			}
			assert(false) : ("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			return 0;
		}

		Arrays.fill(array1, 0, Tools.min(width, maxWidth)+1, big);
		Arrays.fill(array2, 0, Tools.min(width, maxWidth)+1, big);
		arrayCurrent=array1;
		arrayPrev=array2;
		
		{
			if(verbose){System.err.println("\nFirst row.");}
			final byte q=query[qloc];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			edits=big;
			int mloc=1+width-(colLimit-rsloc);
//			assert(false) : width+", "+maxEdits+", "+colLimit+", "+rsloc;
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colLimit-1; col>=colStart; mloc++, col--){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int score=(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; prev=0; score="+score+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			row++; qloc--; rsloc--;
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}
		
		for(row=1; row<len; row++, qloc--, rsloc--){
//			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayCurrent));}
			arrayTemp=arrayCurrent;
			arrayCurrent=arrayPrev;
			arrayPrev=arrayTemp;
			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayPrev)+", current="+Arrays.toString(arrayCurrent));}
			final byte q=query[qloc];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			Arrays.fill(arrayCurrent, big);
			edits=big;
			int mloc=1+width-(colLimit-rsloc);
			boolean forceDiag=(row==len-1);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colLimit-1; col>=colStart; mloc++, col--){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int scoreUp=arrayPrev[mloc+1]+1;
				final int scoreDiag=arrayPrev[mloc]+(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				final int scoreLeft=arrayCurrent[mloc-1]+1;
				final int score=(forceDiag || col==0) ? scoreDiag : Tools.min(scoreUp, scoreDiag, scoreLeft);
				if(verbose){System.err.println("prev=min(s["+(mloc-1)+"]="+arrayCurrent[mloc-1]+", p["+(mloc)+"]="+arrayPrev[mloc]+", p["+(mloc+1)+"]="+arrayPrev[mloc+1]+")");}
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; up="+scoreUp+"; diag="+scoreDiag+"; left="+scoreLeft+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			if(edits>maxEdits){row++; break;}
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}

		lastRow=row-1;
		lastEdits=edits;
		lastOffset=lastOffset(arrayCurrent, halfWidth);
		lastQueryLoc=qloc+1;
		lastRefLoc=rsloc+halfWidth+lastOffset+1;
		while(lastRefLoc<0 || lastQueryLoc<0){lastRefLoc++; lastQueryLoc++;}
		if(verbose){
			System.out.println("\nFinal state: arrayCurrent="+Arrays.toString(arrayCurrent)+"\nlastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
					", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+", qloc="+qloc+", rsloc="+rsloc+"\n");
		}
		return edits;
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	@Override
	public int alignReverseRC(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact){
		assert(big>maxEdits);
		if(verbose){System.err.println("alignReverseRC("+new String(query)+", "+new String(ref)+", "+qstart+", "+rstart+", "+maxEdits+")");}
		if(query.length-qstart>rstart+1){
			int x=alignForwardRC(ref, query, rstart, qstart, maxEdits, exact);
			int temp=lastQueryLoc;
			lastQueryLoc=lastRefLoc;
			lastRefLoc=temp;
			if(verbose){
				System.out.println("Reversed.");
				System.out.println("Final state: lastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
						", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+"\n");
			}
			return x;
		}
		int edits=0, row=0;
		lastRow=-1;
		lastEdits=0;
		lastOffset=0;
		
		final int width=Tools.min(maxWidth, (maxEdits*2)+1, Tools.max(query.length, ref.length)*2+2)|1;
		final int halfWidth=width/2;
		final boolean inexact=!exact;
		
		int qloc=qstart;
		int rsloc=rstart-halfWidth;
		final int xlines=query.length-qstart;
		final int ylines=rstart+1;
		final int len=Tools.min(xlines, ylines);
		if(verbose){System.err.println("xlines="+xlines+", ylines="+ylines+", len="+len);}
		if(len<1){
			if(false){
				throw new RuntimeException("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			}
			assert(false) : ("No overlap: qstart="+qstart+", rstart="+rstart+", qlen="+query.length+", rlen="+ref.length);
			return 0;
		}
		
		Arrays.fill(array1, 0, Tools.min(width, maxWidth)+1, big);
		Arrays.fill(array2, 0, Tools.min(width, maxWidth)+1, big);
		arrayCurrent=array1;
		arrayPrev=array2;
		
		{
			if(verbose){System.err.println("\nFirst row.");}
			final byte q=AminoAcid.baseToComplementExtended[query[qloc]];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			edits=big;
			int mloc=1+width-(colLimit-rsloc);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colLimit-1; col>=colStart; mloc++, col--){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int score=(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; prev=0; score="+score+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			row++; qloc++; rsloc--;
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}
		
		for(row=1; row<len; row++, qloc++, rsloc--){
//			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayCurrent));}
			arrayTemp=arrayCurrent;
			arrayCurrent=arrayPrev;
			arrayPrev=arrayTemp;
			if(verbose){System.err.println("\nNew row, prev="+Arrays.toString(arrayPrev)+", current="+Arrays.toString(arrayCurrent));}
			final byte q=AminoAcid.baseToComplementExtended[query[qloc]];
			final int colStart=Tools.max(0, rsloc);
			final int colLimit=Tools.min(rsloc+width, ref.length);
			Arrays.fill(arrayCurrent, big);
			edits=big;
			int mloc=1+width-(colLimit-rsloc);
			boolean forceDiag=(row==len-1);
			if(verbose){System.err.println("q="+(char)q+", qloc="+qloc+", rsloc="+rsloc+", colStart="+colStart+", colLimit="+colLimit+", mloc="+mloc);}
			for(int col=colLimit-1; col>=colStart; mloc++, col--){
				if(verbose){System.err.println("col="+col+", mloc="+mloc);}
				final byte r=ref[col];
				final int scoreUp=arrayPrev[mloc+1]+1;
				final int scoreDiag=arrayPrev[mloc]+(q==r || (inexact && (!AminoAcid.isFullyDefined(q) || !AminoAcid.isFullyDefined(r))) ? 0 : 1);
				final int scoreLeft=arrayCurrent[mloc-1]+1;
				final int score=(forceDiag || col==0) ? scoreDiag : Tools.min(scoreUp, scoreDiag, scoreLeft);
				if(verbose){System.err.println("prev=min(s["+(mloc-1)+"]="+arrayCurrent[mloc-1]+", p["+(mloc)+"]="+arrayPrev[mloc]+", p["+(mloc+1)+"]="+arrayPrev[mloc+1]+")");}
				arrayCurrent[mloc]=score;
				edits=Tools.min(edits, score);
				if(verbose){System.err.println("Comparing "+(char)q+" to "+(char)r+"; up="+scoreUp+"; diag="+scoreDiag+"; left="+scoreLeft+"; scores = "+Arrays.toString(arrayCurrent));}
			}
			if(edits>maxEdits){row++; break;}
		}
		if(penalizeOffCenter){edits=penalizeOffCenter(arrayCurrent, halfWidth);}
		
		lastRow=row-1;
		lastEdits=edits;
		lastOffset=lastOffset(arrayCurrent, halfWidth);
		lastQueryLoc=qloc-1;
		lastRefLoc=rsloc+halfWidth+lastOffset+1;
		while(lastRefLoc<0 || lastQueryLoc>=query.length){lastRefLoc++; lastQueryLoc--;}
		if(verbose){
			System.out.println("\nFinal state: arrayCurrent="+Arrays.toString(arrayCurrent)+"\nlastRow="+lastRow+", lastEdits="+lastEdits+", lastOffset="+lastOffset+
					", lastQueryLoc="+lastQueryLoc+", lastRefLoc="+lastRefLoc+(query.length<30 ? ", query="+new String(query)+", ref="+new String(ref) : "")+"\n");
		}
		return edits;
	}
	
	/** First working array for dynamic programming matrix calculations */
	private final int[] array1;
	/** Second working array for dynamic programming matrix calculations */
	private final int[] array2;
	private int[] arrayCurrent, arrayPrev, arrayTemp;
	
}
