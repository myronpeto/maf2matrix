package org.panda.resource.tcga;

import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Reads and serves a TCGA mutation file.
 *
 * @author Ozgun Babur
 */
public class MutationReader
{
	private Map<String, Map<String, List<MutTuple>>> mutMap;

	private Set<String> sampleSet;

	public MutationReader(String filename) throws IOException
	{
		this(filename, null);
	}

	public MutationReader(String filename, String... mutTypes) throws IOException
	{
		this.mutMap = new LinkedHashMap<>();
		this.sampleSet = new HashSet<>();
		if (filename != null) load(filename, mutTypes == null || mutTypes.length == 0 ? null : new HashSet<>(Arrays.asList(mutTypes)));
	}

	public void load(String filename, Set<String> mutTypes) throws IOException
	{
		int typeInd = -1;
		int sampleInd = -1;
		int protChInd = -1;

		Optional<String> opt = Files.lines(Paths.get(filename)).filter(l -> !l.startsWith("#"))
			.filter(l -> l.startsWith("Hugo_Symbol")).findFirst();

		if (opt.isPresent())
		{
			String[] header = opt.get().split("\t");
			typeInd = indexOf(header, "Variant_Classification");
			sampleInd = indexOf(header, "Tumor_Sample_Barcode");
			protChInd = indexOf(header, "Protein_Change");
			if (protChInd < 0) protChInd = indexOf(header, "amino_acid_change_WU");
			if (protChInd < 0) protChInd = indexOf(header, "AAChange");
			if (protChInd < 0) protChInd = indexOf(header, "amino_acid_change");
			if (protChInd < 0) protChInd = indexOf(header, "HGVSp_Short");

			if (protChInd < 0)
			{
				System.out.println("No protein change in file " + filename);
//				return;
			}

		}

		processLines(filename, typeInd, sampleInd, protChInd, mutTypes);
	}

	private void processLines(String filename, int typeInd, int sampleInd, int protChInd,
		Set<String> mutTypes) throws IOException
	{
		System.out.println("filename is: " + filename);
		Files.lines(Paths.get(filename)).filter(l -> !l.startsWith("#")).filter(l -> !l.startsWith("Hugo_Symbol"))
			.map(l -> l.split("\t"))
			.filter(t -> !t[0].isEmpty() && !t[0].equals("."))
			.filter(t -> mutTypes == null || mutTypes.contains(t[typeInd]))
			.forEach(token ->
		{
			String id = token[0];
			String sample = token[sampleInd];
			sampleSet.add(sample);

			String type = token[typeInd];

			String protCh = protChInd < 0 || token.length <= protChInd ? "" : token[protChInd];
			if (protCh.startsWith("p.")) protCh = protCh.substring(2);
			else if (protCh.equals(".") || protCh.equals("NULL")) protCh = "";

			MutTuple mut = new MutTuple(type, protCh);

			if (!mutMap.containsKey(id)) mutMap.put(id, new HashMap<>());
			if (!mutMap.get(id).containsKey(sample)) mutMap.get(id).put(sample, new ArrayList<>());
			mutMap.get(id).get(sample).add(mut);
		});
	}

	private boolean multiCenter (String val)
	{
		return Arrays.stream(val.split("\\|")).distinct().count() > 1;
	}

	public Set<String> getSamples()
	{
		return sampleSet;
	}

	public Set<String> getGenes()
	{
		return mutMap.keySet();
	}

	private int indexOf(String[] array, String val)
	{
		for (int i = 0; i < array.length; i++)
		{
			if (array[i].equals(val)) return i;
		}
		return -1;
	}

	/**
	 * All samples have to be in this dataset. This method does not support "no data" conditions.
	 */
	public boolean[] getGeneAlterationArray(String id, String[] samples)
	{
		if (mutMap.containsKey(id))
		{
			boolean[] b = new boolean[samples.length];
			Arrays.fill(b, false);
			for (int i = 0; i < samples.length; i++)
			{
				if (!sampleSet.contains(samples[i]))
					throw new IllegalArgumentException("Sample " + samples[i] + " does not have mutation data.");

				if (mutMap.get(id).containsKey(samples[i])) b[i] = true;
			}
			return b;
		}
		return null;
	}

	/**
	 * @return Array of mutation tuples list. Returns null if id is not recognized. If a sample is not recognized, the
	 * array contains null. An empty list as array element means no mutations in that sample. Do not modify the returned
	 * lists. They are original.
	 */
	public List<MutTuple>[] getMutations(String id, String[] samples)
	{
		if (!mutMap.containsKey(id)) return null;

		List<MutTuple>[] list = new List[samples.length];

		Map<String, List<MutTuple>> stm = mutMap.get(id);

		for (int i = 0; i < samples.length; i++)
		{
			list[i] = stm.get(samples[i]);
			if (list[i] == null && sampleSet.contains(samples[i])) list[i] = Collections.emptyList();
		}

		return list;
	}


	public void writeAsAlterationMatrix(String outFile) throws IOException
	{
		List<String> samples = new ArrayList<>(getSamples());
		Collections.sort(samples);
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

		for (String sample : samples)
		{
			writer.write("\t" + sample);
		}

		for (String gene : mutMap.keySet())
		{
			writer.write("\n" + gene);
			Map<String, List<MutTuple>> map = mutMap.get(gene);

			for (String sample : samples)
			{
				writer.write("\t" + (map.containsKey(sample) ? "1" : "0"));
			}
		}

		writer.close();
	}

	private void printRecurrenceCounts()
	{
		int totalMut = 0;
		int delMut = 0;
		Map<String, Map<String, Integer>> cnt = new HashMap<>();
		for (String gene : mutMap.keySet())
		{
			for (String sample : mutMap.get(gene).keySet())
			{
				List<MutTuple> list = mutMap.get(gene).get(sample);
				for (MutTuple mut : list)
				{
					totalMut++;
					if (mut.value.contains("*") || mut.value.contains("fs")) delMut++;
					if (!cnt.containsKey(gene)) cnt.put(gene, new HashMap<>());
					if (cnt.get(gene).containsKey(mut.value)) cnt.get(gene).put(mut.value, cnt.get(gene).get(mut.value) + 1);
					else cnt.get(gene).put(mut.value, 1);
				}
			}
		}

		System.out.println("Global ratio of deleterious mutations = " + (delMut / (double) totalMut));

		final Map<String, Integer> best = new HashMap<>();
		for (String gene : cnt.keySet())
		{
			best.put(gene, Summary.max(cnt.get(gene).values()));
		}

		List<String> genes = new ArrayList<>(cnt.keySet());
		Collections.sort(genes, (o1, o2) -> best.get(o2).compareTo(best.get(o1)));

		Map<String, Double> dRat = getRatiosOfDeleteriousMutations();

		Histogram h = new Histogram(0.1);
		h.setBorderAtZero(true);

		for (String gene : genes)
		{
			Integer c = best.get(gene);
			if (c == 2) break;
			System.out.println(c + "\t" + gene + "\t" + dRat.get(gene));
			h.count(dRat.get(gene));
		}
		h.print();
	}

	public Map<String, Integer> getHighestRecurrenceCounts()
	{
		Map<String, Map<String, Integer>> cnt = new HashMap<>();
		for (String gene : mutMap.keySet())
		{
			for (String sample : mutMap.get(gene).keySet())
			{
				List<MutTuple> mutTuples = mutMap.get(gene).get(sample);
				for (MutTuple mut : mutTuples)
				{
					if (!cnt.containsKey(gene)) cnt.put(gene, new HashMap<>());
					if (cnt.get(gene).containsKey(mut.value)) cnt.get(gene).put(mut.value, cnt.get(gene).get(mut.value) + 1);
					else cnt.get(gene).put(mut.value, 1);
				}
			}
		}

		Map<String, Integer> highest = new HashMap<>();
		for (String gene : cnt.keySet())
		{
			highest.put(gene, Summary.max(cnt.get(gene).values()));
		}
		return highest;
	}

	public Map<String, Double> getRatiosOfDeleteriousMutations()
	{
		Map<String, Double> rat = new HashMap<>();
		for (String gene : mutMap.keySet())
		{
			int total = 0;
			int del = 0;
			for (String sample : mutMap.get(gene).keySet())
			{
				List<MutTuple> list = mutMap.get(gene).get(sample);
				for (MutTuple mut : list)
				{
					total++;
					if (mut.isDeleterious()) del++;
				}
			}
			double r = del / (double) total;
			rat.put(gene, r);
		}
		return rat;
	}

	public double getOverallDelMutRatio()
	{
		int total = 0;
		int del = 0;

		for (String gene : mutMap.keySet())
		{
			for (String sample : mutMap.get(gene).keySet())
			{
				List<MutTuple> list = mutMap.get(gene).get(sample);
				for (MutTuple mut : list)
				{
					total++;
					if (mut.isDeleterious()) del++;
				}
			}
		}
		double r = del / (double) total;
		return r;
	}

	public Map<String, Integer> getMutatedSampleCounts()
	{
		Map<String, Integer> cnt = new HashMap<>();
		for (String gene : mutMap.keySet())
		{
			cnt.put(gene, mutMap.get(gene).keySet().size());
		}
		return cnt;
	}

	public static void main(String[] args) throws IOException
	{
//		String dir = "/home/ozgun/Documents/TCGA/PanCan/";
//		MutationReader reader = new MutationReader(dir + "tcga_pancancer_082115.vep.filter_whitelisted.maf");
//		reader.writeAsAlterationMatrix(dir + "DataMatrix.txt");

		MutationReader reader = new MutationReader("/home/exacloud/lustre1/users/peto/GDAN/mutation.maf");
//		reader.printRecurrenceCounts();
	}
}
