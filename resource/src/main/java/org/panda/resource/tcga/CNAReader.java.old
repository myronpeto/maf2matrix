package org.panda.resource.tcga;

import org.panda.utility.ArrayUtil;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.TTest;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Reads and serves a TCGA Gistic Copy Number file.
 *
 * @author Ozgun Babur
 */
public class CNAReader
{
	private String filename;

	private Map<String, Map<String, Integer>> data;

	boolean reduce;

	int threshold;

	public static final int NO_DATA = -Integer.MAX_VALUE;

	public CNAReader(String filename) throws FileNotFoundException
	{
		this(filename, null);
	}

	public CNAReader(String filename, Set<String> genes) throws FileNotFoundException
	{
		this(filename, genes, true, 2);
	}

	public CNAReader(String filename, boolean reduce, int threshold) throws FileNotFoundException
	{
		this(filename, null, reduce, threshold);
	}

	public CNAReader(String filename, Set<String> genes, boolean reduce, int threshold) throws FileNotFoundException
	{
		this.filename = filename;
		this.reduce = reduce;

		if (reduce && threshold < 1)
			throw new IllegalArgumentException("Threshold has to be positive integer");

		this.threshold = threshold;
		this.data = new HashMap<String, Map<String, Integer>>();
		load(genes);
	}

	private void load(Set<String> genes) throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(filename));
		String line = sc.nextLine();
		while (line.startsWith("#")) line = sc.nextLine();

		String[] header = line.split("\t");

		int ss = 0;
		while (!header[ss].startsWith("TCGA")) ss++;

		for (int i = ss; i < header.length; i++)
		{
			header[i] = header[i].substring(0, 12);
		}

		while (sc.hasNextLine())
		{
			line = sc.nextLine();
			String id = line.substring(0, line.indexOf("\t"));
			if (id.contains("|")) id = id.substring(0, id.indexOf("|"));

			if (genes != null && !genes.contains(id)) continue;

			String[] token = line.split("\t");


			for (int i = ss; i < header.length; i++)
			{
				Integer val = Integer.parseInt(token[i]);

//				if (val != 0)
				{
					if (!data.containsKey(id)) data.put(id, new HashMap<String, Integer>());
					data.get(id).put(header[i], val);
				}
			}

			if (genes != null && genes.size() == data.size()) break;
		}
	}

	public Set<String> getSamples()
	{
		Set<String> samples = new HashSet<String>();
		for (Map<String, Integer> map : data.values())
		{
			samples.addAll(map.keySet());
		}
		return samples;
	}

	public Set<String> getGenes()
	{
		return data.keySet();
	}

	public int[] getGeneAlterationArray(String id, String[] samples)
	{
		if (data.containsKey(id))
		{
			int[] b = new int[samples.length];
			Arrays.fill(b, 0);
			for (int i = 0; i < samples.length; i++)
			{
				if (data.get(id).containsKey(samples[i]))
				{
					Integer val = data.get(id).get(samples[i]);
					if (reduce) b[i] = val >= threshold ? 1 : val <= -threshold ? -1 : 0;
					else b[i] = val;
				}
				else
				{
					b[i] = NO_DATA;
				}
			}
			return b;
		}
		return null;
	}

	private boolean[] getAmplified(int[] alterations)
	{
		return getOneSidedChange(alterations, true);
	}

	private boolean[] getDeleted(int[] alterations)
	{
		return getOneSidedChange(alterations, false);
	}

	private boolean[] getOneSidedChange(int[] alterations, boolean amplified)
	{
		boolean[] b = new boolean[alterations.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = (amplified && alterations[i] > 0) || (!amplified && alterations[i] < 0);
		}
		return b;
	}

	private boolean[] getNoChange(int[] alterations)
	{
		boolean[] b = new boolean[alterations.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = alterations[i] == 0;
		}
		return b;
	}

	public int[] getExpVerifiedCNA(String id, String[] samples, double[] exp, double pvalThr)
	{
		assert exp.length == samples.length;
		int[] alterations = getGeneAlterationArray(id, samples);

		boolean[] noChange = getNoChange(alterations);
		double[] noChVals = ArrayUtil.subset(exp, noChange);
		if (noChVals.length == 0) return null;

		boolean[] amplified = getAmplified(alterations);
		boolean[] deleted = getDeleted(alterations);

		amplified = getVerified(amplified, noChVals, true, exp, pvalThr);
		deleted = getVerified(deleted, noChVals, false, exp, pvalThr);

		if (amplified == null && deleted == null) return null;
		if (amplified == null) return keepSelected(alterations, deleted);
		if (deleted == null) return keepSelected(alterations, amplified);

		ArrayUtil.ORWith(amplified, deleted);
		return keepSelected(alterations, amplified);
	}

	public boolean[] getVerified(boolean[] change, double[] noChVals, boolean amplified,
		double[] exp, double pvalThr)
	{
		double[] chVals = ArrayUtil.subset(exp, change);
		if (chVals.length == 0) return null;

		double ch = Summary.calcChangeOfMean(noChVals, chVals);
		if ((amplified && ch < 0) || (!amplified && ch > 0))
		{
			return null;
		}

		double p = TTest.getPValOfMeanDifference(noChVals, chVals);
		if (Double.isNaN(p) || p > pvalThr) return null;

		double val = Summary.getIntersectionPoint(chVals, noChVals);

		return select(change, exp, val, amplified);
	}

	private static boolean[] select(boolean[] considerLoc, double[] exp, double thr, boolean greaterThan)
	{
		boolean[] select = new boolean[exp.length];

		for (int i = 0; i < select.length; i++)
		{
			select[i] = considerLoc[i] && (greaterThan ? exp[i] > thr : exp[i] < thr);
		}
		return select;
	}

	private int[] keepSelected(int[] cna, boolean[] verified)
	{
		int[] arr = new int[cna.length];
		for (int i = 0; i < arr.length; i++)
		{
			arr[i] = verified[i] ? cna[i] : 0;
		}
		return arr;
	}

	public double getCNAPval(String id, String[] samples, boolean amplified, double[] exp)
	{
		int[] alterations = getGeneAlterationArray(id, samples);

		boolean[] noChange = getNoChange(alterations);
		double[] noChVals = ArrayUtil.subset(exp, noChange);
		if (noChVals.length == 0) return 1;

		boolean[] changed = amplified ? getAmplified(alterations) : getDeleted(alterations);
		double[] chVals = ArrayUtil.subset(exp, changed);

		if (chVals.length == 0) return 1;

		double ch = Summary.calcChangeOfMean(noChVals, chVals);
		if ((amplified && ch < 0) || (!amplified && ch > 0))
		{
			return 1;
		}

		double p = TTest.getPValOfMeanDifference(noChVals, chVals);
		if (Double.isNaN(p)) return 1;

		return p;
	}

	public static void main(String[] args) throws FileNotFoundException
	{
		CNAReader reader = new CNAReader("/home/ozgun/Documents/TCGA/UVM/all_thresholded.by_genes.txt");
		System.out.println(reader.getSamples().size());
	}
}
