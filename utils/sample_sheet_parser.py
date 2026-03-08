from abc import ABC, abstractmethod
from argparse import ArgumentParser, Namespace
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, LiteralString
import numpy as np
import pandas as pd


# ! To get rid of spelling issues,
# ! only lower case names and hyphen-minus (-) are allowed.
# ! Also, headers only contain upper case letters and underscore (_) to distinguish.
PIPELINE_CODING: dict[str, dict[str, str]] = {
    'TRIMMER': {
        'no-trim': '0', 'fastp': '1', 'trim-galore': '2'
    },
    'ALIGNER': {
        'bismark-bowtie2': '0', 'bismark-hisat2': '1',
        'hisat2': '2', 'bowtie2': '3',
        'bwa-meth': '4', 'bwa-meme': '5', 'dnmtools': '6',
        'biscuit': '7', 'bsmapz': '8', 'batmeth2': '9',
        'bsbolt': 'A', 'abismal': 'B', 'hisat-3n': 'C',
        'bitmapperbs': 'D', 'whisper': 'E', 'pufferfish': 'F',
        'strobealign': 'G', 'aryana': 'H', 'gem3': 'I',
        'segemehl': 'J', 'last': 'K'
    },
    'DEDUPER': {
        'no-dedup': '0', 'bismark-dedup': '1', 'gatk-dedup': '2',  # .../gatk/gatk/... will cause ambiguity.
        'samtools': '3', 'sambamba': '4', 'nubeam-dedup': '5',
        'bio-seq-zip': '6', 'trie-dedup': '7', 'ngs-reads-treatment': '8',
        'dupsifter': '9', 'gencore': 'A', 'minirmd': 'B', 'gatk-optical': 'C'
    },
    'CALIBRATOR': {
        'no-pre-calibrate': '@', 'no-calibrate': '0', 'gatk-cali': '1', 'gatk': '1', 'care': '2',
        'lighter': '3', 'reckoner': '4'
    },
    'COUNTER': {
        'bismark': '0', 'biscuit': '1', 'methyldackel': '2',
        'dnmtools': '3', 'astair': '4', 'rastair': '5',
        'bsgenova': '6', 'fame': '7', 'deepvariant': '8',
        'haplotypecaller': '9', 'methylationtypecaller': 'A',
        'epiallele': 'B'
    },  # row for FAME: BASENAME,TRIMMER,FASTQCER,,,,fame,False
    'FASTQCER': {
        'fastqc': '0', 'falco': '1'  # actually FASTQCER and BAMSTATIST are not coded,
    },  # just put it here to make sure the right one included.
    'BAMSTATIST': {'1': '1'}
}

rev_pipeline_coding: dict[str, dict[str, str]] = {k: {v: k for k, v in v.items()} for k, v in PIPELINE_CODING.items()}

ENCODING_ORDER: list[str] = ['TRIMMER', 'DEDUPER', 'CALIBRATOR',
                             'ALIGNER', 'DEDUPER', 'CALIBRATOR', 'COUNTER',]

DEDUPER_BEFORE_ALIGNMENT: set[str] = {'nubeam-dedup', 'bio-seq-zip', 'trie-dedup', 'ngs-reads-treatment', 'minirmd'}
CALIBRATOR_BEFORE_ALIGNMENT: set[str] = {'bfc', 'care', 'lighter', 'reckoner', 'sequencerr'}


@dataclass
class PipelineContext:
    fname: str
    parent_dir: Path
    trimmed: bool = False
    aligned: bool = False
    deduped: bool = False
    calibrated: bool = False
    counted: bool = False
    trace: list[str] = field(default_factory=list)  # type: ignore

    def update(self,
               parent_dir: Path | None = None,
               trimmed: bool | None = None,
               aligned: bool | None = None,
               deduped: bool | None = None,
               calibrated: bool | None = None,
               counted: bool | None = None) -> 'PipelineContext':
        if sum(1 for attr in [trimmed, aligned, deduped, calibrated, counted] if attr is not None) > 1:
            raise ValueError('Only one operation at a time.')

        new_trace = self.trace.copy()
        if trimmed is True and not self.trimmed:
            new_trace.append('trimmed')
        if aligned is True and not self.aligned:
            new_trace.append('aligned')
        if deduped is True and not self.deduped:
            new_trace.append('deduped')
        if calibrated is True and not self.calibrated:
            new_trace.append('calibrated')
        if counted is True and not self.counted:
            new_trace.append('counted')

        return PipelineContext(fname=self.fname,
                               parent_dir=parent_dir if parent_dir is not None else self.parent_dir,
                               trimmed=trimmed if trimmed is not None else self.trimmed,
                               aligned=aligned if aligned is not None else self.aligned,
                               deduped=deduped if deduped is not None else self.deduped,
                               calibrated=calibrated if calibrated is not None else self.calibrated,
                               counted=counted if counted is not None else self.counted,
                               trace=new_trace)

    def done_before(self, module1: str, module2: str) -> bool:
        for m in (module1, module2):
            if m not in self.trace:
                raise ValueError(f'{m} not done yet.')

        if self.trace.index(module1) < self.trace.index(module2):
            return True
        else:
            return False


@dataclass
class EmptyModule:
    module: str

    def __bool__(self):
        return False


class PipelineModule(ABC):
    @property
    @abstractmethod
    def tool(self) -> str:
        pass

    def __eq__(self, other):  # type: ignore
        if isinstance(other, self.__class__):
            return self.tool == other.tool
        else:
            return False

    def __str__(self):
        return f'{self.tool} as {self.__class__.__name__}'

    def __hash__(self) -> int:
        return hash((type(self), self.tool))

    @abstractmethod
    def files(self, context: PipelineContext) -> tuple[list[Path], PipelineContext]:
        pass


class Trimmer(PipelineModule):
    # ! YEAH, the file names spawn here.
    EXT_MAPPING: dict[str, list[str]] = {
        'fastp': ['.R1.fq.gz', '.R2.fq.gz',
                  '.fastp.json', '.fastp.html'],
        'trim-galore': ['.R1.fq.gz', '.R2.fq.gz',
                        '.R1.report', '.R2.report'],
        'no-trim': ['.R1.fq.gz', '.R2.fq.gz']
    }

    def __init__(self, trimmer: str):
        self._tool: str = trimmer

    @property
    def tool(self) -> str:
        return self._tool

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.tool}) '
                f'-> files: {self.EXT_MAPPING[self.tool]}')

    def __str__(self):
        if self.tool == 'no-trim':
            return 'Skipped: no Trimmer'
        else:
            return super().__str__()

    def files(self, context: PipelineContext) -> tuple[list[Path], PipelineContext]:
        return ([context.parent_dir / self.tool / f'{context.fname}{ext}'
                 for ext in self.EXT_MAPPING[self.tool]],
                context.update(parent_dir=context.parent_dir / self.tool,
                               trimmed=True))


class Fastqcer(PipelineModule):
    EXT_MAPPING: dict[str, list[str]] = {
        'fastqc': ['.R1.fastqc.html', '.R2.fastqc.html'],
        'falco': ['.R1.falco.html', '.R2.falco.html']}

    def __init__(self, fastqcer: str):
        self._tool: str = fastqcer

    @property
    def tool(self) -> str:
        return self._tool

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.tool})'
                f' -> files: {self.EXT_MAPPING[self.tool]}')

    def files(self, context: PipelineContext) -> tuple[list[Path], PipelineContext]:
        return ([context.parent_dir / f'{context.fname}{ext}'
                 for ext in self.EXT_MAPPING[self.tool]],
                context)


class Aligner(PipelineModule):
    def __init__(self, aligner: str):
        self._tool: str = aligner

    @property
    def tool(self) -> str:
        return self._tool

    def __repr__(self):
        return f'{self.__class__.__name__}({self.tool}) -> files: [.bam", ".bam.bai"]'

    def files(self, context: PipelineContext) -> tuple[list[Path], PipelineContext]:
        current_dir: Path
        files: list[Path]

        if context.deduped and not context.calibrated:
            current_dir = context.parent_dir / 'no-pre-calibrate' / self.tool
        elif not context.deduped and not context.calibrated:
            current_dir = context.parent_dir / 'no-pre-dedup' / 'no-pre-calibrate' / self.tool
        else:
            current_dir = context.parent_dir / self.tool

        files = [current_dir / f'{context.fname}{ext}'
                 for ext in ('.bam', '.bam.bai')]
        if self.tool in ('msuite2-bowtie2',
                         'msuite2-hisat2'):
            files += [current_dir / f
                      for f in (f'{context.fname}.msuite2.bedgraph.gz',
                                'Msuite2.report/index.html')]
        return files, context.update(parent_dir=current_dir, aligned=True)


class DeDuper(PipelineModule):
    EXT_MAPPING = {
        'no-dedup': ['.bam', '.bam.bai'],
        'bismark-dedup': ['.bam', '.bam.bai'],
        'gatk-dedup': ['.bam', '.bam.bai', '.dedup.metrics'],
        'samtools': ['.bam', '.bam.bai', '.dup.stats'],
        'sambamba': ['.bam', '.bam.bai'],
        'nubeam-dedup': ['.R1.fq.gz', '.R2.fq.gz'],
        'bio-seq-zip': ['.R1.fq.gz', '.R2.fq.gz', '.bioseqzip.report'],
        'trie-dedup': ['.R1.fq.gz', '.R2.fq.gz'],
        'ngs-reads-treatment': ['.R1.fq.gz', '.R2.fq.gz'],
        'dupsifter': ['.bam', '.bam.bai', '.dupsifter.stat'],
        'gencore': ['.bam', '.bam.bai', '.gencore.json', 'gencore.html'],
        'minirmd': ['.R1.fq.gz', '.R2.fq.gz'],
        'gatk-optical': ['.bam', '.bam.bai', '.dedup.metrics']
    }

    def __init__(self, deduper: str):
        self._tool: str = deduper

    @property
    def tool(self) -> str:
        return self._tool

    def __repr__(self):
        return f'{self.__class__.__name__}({self.tool}) -> files: {self.EXT_MAPPING[self.tool]}'

    def __str__(self):
        if self.tool == 'no-dedup':
            return 'Skipped: no DeDuper'
        else:
            return super().__str__()

    def files(self, context: PipelineContext) -> tuple[list[Path], PipelineContext]:
        return ([context.parent_dir / self.tool / f'{context.fname}{ext}'
                 for ext in self.EXT_MAPPING[self.tool]],
                context.update(parent_dir=context.parent_dir / self.tool,
                               deduped=True))


class Calibrator(PipelineModule):
    EXT_MAPPING = {
        'no-calibrate': ['.bam', '.bam.bai'],
        'care': ['.R1.fq.gz', '.R2.fq.gz'],
        'lighter': ['.R1.fq.gz', '.R2.fq.gz'],
        'gatk-cali': ['.bam', '.bam.bai',
                      '.before.table', '.after.table',
                      '.AnalyzeCovariates.pdf'],
        'reckoner': ['.R1.fq.gz', '.R2.fq.gz'],
        'sequencerr': ['.R1.fq.gz', '.R2.fq.gz']
    }

    def __init__(self, calibrator: str):
        self._tool: str = calibrator

    @property
    def tool(self) -> str:
        return self._tool

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.tool})'
                f' -> files: {self.EXT_MAPPING[self.tool]}')

    def __str__(self):
        if self.tool == 'no-calibrate':
            return 'Skipped: no Calibrator'
        else:
            return super().__str__()

    def files(self, context: PipelineContext) -> tuple[list[Path], PipelineContext]:
        current_dir: Path

        if context.deduped and context.aligned and context.done_before('deduped',
                                                                       'aligned'):
            current_dir = context.parent_dir / 'no-dedup' / self.tool

        elif not context.deduped and not context.aligned:
            current_dir = context.parent_dir / 'no-pre-dedup' / self.tool
        else:
            current_dir = context.parent_dir / self.tool

        return ([current_dir / f'{context.fname}{ext}'
                 for ext in self.EXT_MAPPING[self.tool]],
                context.update(parent_dir=current_dir, calibrated=True))


class Counter(PipelineModule):
    EXT_MAPPING = {
        'bismark': ['.bedgraph.gz', '.ucsc_bedgraph.gz',
                    '.CHG_OB.txt.gz', '.CHG_OT.txt.gz',
                    '.CHH_OB.txt.gz', '.CHH_OT.txt.gz',
                    '.CpG_OB.txt.gz', '.CpG_OT.txt.gz',
                    '.bismark.cov.gz', '.splitting_report',
                    '.c2c.cov.gz', '.c2c.report.gz', '.c2c.summary'],
        'astair': ['.astair.mods.gz', '.astair.stats'],
        'rastair': ['.rastair.mods.gz'],
        'dnmtools': ['.dnmtools.gz', '.dnmtools.stats'],
        'methyldackel': ['.bedgraph.gz', '.merged.bedgraph.gz'],
        'biscuit': ['.epibed.gz'],
        'bsgenova': ['.ATCGmap.gz', '.CGmap.gz', '.bed.gz'],
        'fame': ['.bedgraph.zst'],
        'deepvariant': ['.deepvariant.vcf.gz', '.deepvariant.gvcf.gz',
                        '.deepvariant.bedgraph.zst'],
        'haplotypecaller': ['.gvcf.gz', '.vcf.gz', '.bedgraph.zst'],
        'methylationtypecaller': ['.vcf.gz', '.bedgraph.zst'],
        'epiallele': ['.report.zst']
    }

    def __init__(self, counter: str):
        self._tool: str = counter

    @property
    def tool(self) -> str:
        return self._tool

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.tool}) '
                f'-> files: {self.EXT_MAPPING[self.tool]}')

    def files(self, context: PipelineContext) -> tuple[list[Path], PipelineContext]:
        current_dir: Path

        if context.aligned:
            match (context.done_before('calibrated',
                                       'aligned'),
                   context.done_before('deduped',
                                       'aligned')):
                case (True, True):
                    current_dir = context.parent_dir / 'no-dedup' / 'no-calibrate' / self.tool
                case (False, True) | (False, False):
                    current_dir = context.parent_dir / self.tool
                case _:
                    current_dir = context.parent_dir / 'no-calibrate' / self.tool
        else:  # FAME
            match (context.deduped, context.calibrated):
                case (True, False):
                    current_dir = context.parent_dir / 'no-pre-calibrate' / self.tool
                case (False, False):
                    current_dir = context.parent_dir / 'no-pre-dedup' / 'no-pre-calibrate' / self.tool
                case _:
                    current_dir = context.parent_dir / self.tool

        return ([current_dir / f'{context.fname}{ext}'
                 for ext in self.EXT_MAPPING[self.tool]],
                context.update(parent_dir=current_dir, counted=True))


class BAMStatist(PipelineModule):
    EXT_TUPLE = ('.samtools.stats.txt', '.samtools.flagstats.txt',
                 '.methydackel_mbias_OT.svg', '.methydackel_mbias_OB.svg'
                                              '.biscuit.CpGRetentionByReadPos.txt',
                 '.biscuit.CpHRetentionByReadPos.txt',
                 '.biscuit.dup_report.txt', '.biscuit.isize_table.txt',
                 '.biscuit.mapq_table.txt', '.biscuit.strand_table.txt',
                 '.biscuit.totalReadConversionRate.txt',
                 '.biscuit.bsstrand.txt', '.biscuit.cinread.txt',
                 '.qualimap.pdf')

    def __init__(self):
        self._tool: str = 'samtools,methydackel,biscuit,bismark,qualimap'
        self.tool_list: list[str] | list[LiteralString] = self._tool.split(',')

    @property
    def tool(self) -> str:
        return self._tool

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.tool})'
                f' -> files: {self.EXT_TUPLE}')

    def files(self, context: PipelineContext) -> tuple[list[Path], PipelineContext]:
        fname: str = context.fname
        parent_dir: Path = context.parent_dir

        files: list[Path]
        files = [
            parent_dir / f'{fname}{ext}'
            for ext in self.EXT_TUPLE
        ] + [
            parent_dir / 'qualimap' / ext
            for ext in (f'{fname}.qualimap.pdf', 'qualimapReport.html')
        ]
        if any(bismark in parent_dir.as_posix()
               for bismark in ('bismark-bowtie2', 'bismark-hisat2')):
            files.append(parent_dir / f'{fname}.bam2nuc.txt')
        return files, context


RegisteredModule = PipelineModule | EmptyModule
TrimmerOption = Trimmer | EmptyModule
FastqcerOption = Fastqcer | EmptyModule
AlignerOption = Aligner | EmptyModule
DeDuperOption = DeDuper | EmptyModule
CalibratorOption = Calibrator | EmptyModule
CounterOption = Counter | EmptyModule
BAMStatistOption = BAMStatist | EmptyModule
CoreModule = Trimmer | Aligner | DeDuper | Calibrator | Counter


class Fq2BedGraphRun:
    def __init__(self, row: pd.Series,   # type: ignore
                 root_dir: Path):
        self.row: pd.Series = row  # type: ignore
        self.BaseName: str = str(row['FNAME'])  # type: ignore
        self.Trimmer: Trimmer = Trimmer(str(row['TRIMMER']))  # type: ignore
        self.Fastqcer: FastqcerOption = (Fastqcer(str(row['FASTQCER']))  # type: ignore
                                         if self._vaild_str(row['FASTQCER'])  # type: ignore
                                         else EmptyModule('Fastqcer'))
        self.Aligner: AlignerOption = (Aligner(str(row['ALIGNER']))  # type: ignore
                                       if self._vaild_str(row['ALIGNER'])  # type: ignore
                                       else EmptyModule('Aligner'))
        self.DeDuper: DeDuperOption = (DeDuper(str(row['DEDUPER']))  # type: ignore
                                       if self._vaild_str(row['DEDUPER'])  # type: ignore
                                       else EmptyModule('DeDuper'))
        self.Calibrator: CalibratorOption = (Calibrator(str(row['CALIBRATOR']))  # type: ignore
                                             if self._vaild_str(row['CALIBRATOR'])  # type: ignore
                                             else EmptyModule('Calibrator'))
        self.Counter: CounterOption = (Counter(str(row['COUNTER']))  # type: ignore
                                       if self._vaild_str(row['COUNTER'])  # type: ignore
                                       else EmptyModule('Counter'))
        self.BAMStatist: BAMStatistOption = (BAMStatist()
                                             if self._vaild_str(row['STATS'])
                                             else EmptyModule('BAMStatist'))
        self.ModuleChain = self._tool_chain()
        self.ParentDir = self._parent_dir(root_dir=root_dir, fname=self.BaseName)
        self.Files = self._files()
        self.Code = self._chain_encode(pipeline_coding=PIPELINE_CODING)

    def __repr__(self):
        return (f'{self.__class__.__name__} for {self.row["FNAME"]}\n'  # type: ignore
                f'Trimmer: {self.Trimmer}\n'
                f'QCReporter: {self.Fastqcer}\n'
                f'Aligner: {self.Aligner}\n'
                f'DeDuper: {self.DeDuper}\n'
                f'Calibrator: {self.Calibrator}\n'
                f'Counter: {self.Counter}\n'
                f'BAMStatist: {self.BAMStatist}\n'
                f'ModuleChain: {self.ModuleChain}\n'
                f'ParentDir: {self.ParentDir}')

    @staticmethod
    def _vaild_str(some_obj: Any):
        """
        To get rid of annoying pyright error message.
        """
        return isinstance(some_obj, str) and len(some_obj) > 0

    @staticmethod
    def _parent_dir(root_dir: Path, fname: str) -> Path:
        return root_dir / fname

    @staticmethod
    def _order_check(module: RegisteredModule,
                     pre_tools: set[str]) -> Literal[-1, 0, 1]:
        check_result = -1
        if isinstance(module, PipelineModule):
            if module.tool in pre_tools:
                check_result = 1
            else:
                check_result = 0
        return check_result

    def _tool_chain(self) -> list[PipelineModule]:
        full_toolchain: list[RegisteredModule]
        tool_chain: list[PipelineModule]

        if isinstance(self.Counter, PipelineModule) and self.Counter.tool == 'fame':
            for _m, _l in [[self.DeDuper, DEDUPER_BEFORE_ALIGNMENT],
                           [self.Calibrator, CALIBRATOR_BEFORE_ALIGNMENT]]:
                if isinstance(_m, EmptyModule):
                    pass
                else:
                    if _m.tool not in _l:  # type: ignore
                        raise ValueError(f'{_m} not in {_l}')
                    else:
                        pass
            if isinstance(self.DeDuper, PipelineModule) or isinstance(self.Calibrator, PipelineModule):
                full_toolchain = [i for i in [self.Trimmer, self.Fastqcer, self.DeDuper,
                                              self.Calibrator, self.Fastqcer, self.Counter]
                                  if isinstance(i, PipelineModule)]
            else:
                full_toolchain = [self.Trimmer, self.Fastqcer, self.Counter]
        else:
            match (self._order_check(module=self.DeDuper,
                                     pre_tools=DEDUPER_BEFORE_ALIGNMENT),
                   self._order_check(module=self.Calibrator,
                                     pre_tools=CALIBRATOR_BEFORE_ALIGNMENT)):
                # ! -1 = None; 0 = post; 1 = pre.
                # ! When -1, treat the left chain as conventional one.
                case (1, -1) | (1, 0):
                    full_toolchain = [self.Trimmer, self.Fastqcer, self.DeDuper, self.Fastqcer,
                                      self.Aligner, self.BAMStatist, self.Calibrator, self.BAMStatist, self.Counter]
                case (-1, 1) | (0, 1):
                    full_toolchain = [self.Trimmer, self.Fastqcer, self.Calibrator, self.Fastqcer,
                                      self.Aligner, self.BAMStatist, self.DeDuper, self.BAMStatist, self.Counter]
                case (1, 1):
                    full_toolchain = [self.Trimmer, self.Fastqcer, self.DeDuper, self.Calibrator,
                                      self.Fastqcer, self.Aligner, self.BAMStatist, self.Counter]
                case _:
                    full_toolchain = [self.Trimmer, self.Fastqcer, self.Aligner, self.BAMStatist,
                                      self.DeDuper, self.Calibrator, self.BAMStatist, self.Counter]

        tool_chain = []
        for module in full_toolchain:
            if isinstance(module, PipelineModule):
                tool_chain.append(module)
                if module.tool in ('msuite2-bowtie2', 'msuite2-hisat2'):
                    break
            else:
                if module.module not in ('Fastqcer', 'BAMStatist'):
                    break
        return tool_chain

    def _files(self) -> list[Path]:
        files: list[Path] = []
        last_non_stat_files = []
        current_context: PipelineContext = PipelineContext(fname=self.BaseName,
                                                           parent_dir=self.ParentDir)
        for tool in self.ModuleChain:
            (new_files,
             current_context) = tool.files(context=current_context)
            if any(isinstance(tool, module) for module in (Fastqcer, BAMStatist)):
                files += new_files
            else:
                last_non_stat_files = new_files

        files += last_non_stat_files

        return files

    def module_context(self, module_type: str) -> Path:
        module_dir: Path | None = None
        current_context: PipelineContext = PipelineContext(fname=self.BaseName,
                                                           parent_dir=self.ParentDir)
        for module in self.ModuleChain:
            (_,
             current_context) = module.files(context=current_context)
            if module.__class__.__name__ == module_type:
                module_dir = current_context.parent_dir / module.tool
                break
        if not module_dir:
            raise ValueError(f'{module_type} not found in the current {self.ModuleChain}.')
        return module_dir

    @staticmethod
    def _code(module: RegisteredModule, coding_dict: dict[str, dict[str, str]]):
        try:
            return coding_dict[module.__class__.__name__.upper()][module.tool]  # type: ignore
        except (KeyError, AttributeError, ValueError):
            return '?'

    def _pipeline_encode(self,
                         pipeline_coding: dict[str, dict[str, str]] = PIPELINE_CODING,
                         encoding_order: list[str] = ENCODING_ORDER) -> str:
        # ! do not code the ModuleChain. Instead, code the final context.
        if isinstance(self.Counter, EmptyModule) and \
                isinstance(self.Aligner, PipelineModule) and \
                self.Aligner.tool in ('msuite2-bowtie2', 'msuite2-hisat2'):
            # which do not need a counter
            counter_dir_parts = (self.module_context(module_type='Aligner')
                                 .relative_to(self.ParentDir)
                                 .parts)
        else:
            counter_dir_parts = (self.module_context(module_type='Counter')
                                     .relative_to(self.ParentDir)
                                     .parts)
        # trimmer, pre_dedup, pre_calibrate, align, dedup, calibrate, counter = relative_counter_dir.parts
        return ''.join(pipeline_coding[module_type].get(dir_part, '?')
                       for dir_part, module_type in zip(counter_dir_parts, encoding_order))

    def _chain_encode(self,
                      pipeline_coding: dict[str, dict[str, str]] = PIPELINE_CODING,
                      encoding_order: list[str] = ENCODING_ORDER) -> str:
        last_non_stat_files: list[Path] = []
        current_context: PipelineContext = PipelineContext(fname=self.BaseName,
                                                           parent_dir=self.ParentDir)
        for tool in self.ModuleChain:
            (new_files,
             current_context) = tool.files(context=current_context)
            if any(isinstance(tool, module) for module in (Fastqcer, BAMStatist)):
                pass
            else:
                last_non_stat_files = new_files

        final_dir_parts = last_non_stat_files[0].relative_to(self.ParentDir).parts

        return ''.join(pipeline_coding[module_type].get(dir_part, '?')
                       for dir_part, module_type in zip(final_dir_parts, encoding_order))


def pipeline_decode(pipeline_code: str,
                    pipeline_decoding: dict[str, dict[str, str]] = PIPELINE_CODING,
                    encoding_order: list[str] = ENCODING_ORDER) -> Path:
    if '?' in pipeline_code and pipeline_code[3] not in ['2', '3']:  # msuit2-bowtie2 and msuit-hisat2
        raise ValueError('Invalid pipeline code.')

    if pipeline_code[3] in ['2', '3']:  # msuit2-bowtie2 and msuit-hisat2
        return Path('/'.join(pipeline_decoding[module_type][code]
                             for module_type, code in zip(encoding_order, pipeline_code[: 4])))
    else:
        return Path('/'.join(pipeline_decoding[module_type][code]
                             for module_type, code in zip(encoding_order, pipeline_code)))


def read_fq2bedgraph_sheet(csv_path: Path,
                           pipeline_coding: dict = PIPELINE_CODING) -> pd.DataFrame:  # type: ignore
    query_string = ('FNAME.notna() and '
                    '(TRIMMER.notna() | TRIMMER.isin(@pipeline_coding["TRIMMER"].keys())) and '
                    '(FASTQCER.isna() | FASTQCER.isin(@pipeline_coding["FASTQCER"].keys())) and '
                    '(ALIGNER.isna() | ALIGNER.isin(@pipeline_coding["ALIGNER"].keys())) and '
                    '(DEDUPER.isna() | DEDUPER.isin(@pipeline_coding["DEDUPER"].keys())) and '
                    '(CALIBRATOR.isna() | CALIBRATOR.isin(@pipeline_coding["CALIBRATOR"].keys())) and '
                    '(COUNTER.isna() | COUNTER.isin(@pipeline_coding["COUNTER"].keys())) and '
                    '(STATS.isna() | STATS.isin([True, False]))')

    input_df: pd.DataFrame
    f2b_df: pd.DataFrame

    input_df = (pd.read_csv(csv_path, comment='#',  # type: ignore
                            dtype={'FNAME': 'str', 'TRIMMER': 'str',
                                   'FASTQCER': 'str', 'ALIGNER': 'str',
                                   'DEDUPER': 'str', 'CALIBRATOR': 'str',
                                   'COUNTER': 'str', 'STATS': 'boolean'}))

    f2b_df = (input_df.dropna(subset=['FNAME', 'TRIMMER'])  # type: ignore
                      .query(query_string, engine='python'))
    # !LEGALITY CHECK HERE BEFORE INITIALIZATION. NoneType IS LEGAL.

    if f2b_df.empty:
        wrong_dict = {j: [i for i in input_df[j].unique()  # type: ignore
                          if i and i not in pipeline_coding[j].keys()]  # type: ignore
                      for j in ['FASTQCER', 'TRIMMER', 'ALIGNER',
                                'DEDUPER', 'CALIBRATOR', 'COUNTER']}
        wrong_message = 'No valid rows found in the sheet. Diagnosis:'
        for category, wrong_ls in wrong_dict.items():
            if wrong_ls:
                wrong_message += (f'\nwrong {category}:\n{wrong_ls}\n'
                                  f'avaliable:\n{pipeline_coding[category].keys()}')  # type: ignore
            else:
                wrong_message += f'{category} fine'

        if wrong_message == 'No valid rows found in the sheet. Diagnosis:':
            wrong_message += ('all fine: FASTQCER, TRIMMER, ALIGNER, DEDUPER, CALIBRATOR, COUNTER'
                              'check TRIMMER or STATS.')

        raise ValueError(wrong_message)

    return f2b_df.replace(np.nan, None)  # type: ignore


def fq2bedgraph_file_ls(csv_path_str: str) -> list[str]:
    df = read_fq2bedgraph_sheet(csv_path=Path(csv_path_str))
    return ([
        file.as_posix()
        for _, row in df.iterrows()  # type: ignore
        for file in Fq2BedGraphRun(row=row, root_dir=Path('result')).Files
    ])


def dev_fq2bedgraph_file_ls(df: pd.DataFrame) -> list[str]:
    return ([
        file.relative_to(Path()).as_posix()
        for _, row in df.iterrows()  # type: ignore
        for file in Fq2BedGraphRun(row=row, root_dir=Path('result')).Files
    ])


def fq2bedgraph_tool_ls(csv_path_str: str) -> dict[str, list[str] | list[LiteralString]]:
    df = read_fq2bedgraph_sheet(csv_path=Path(csv_path_str))
    tool_chain_set: set[PipelineModule] = {
        tool
        for _, row in df.iterrows()  # type: ignore
        for tool in Fq2BedGraphRun(row=row, root_dir=Path('result')).ModuleChain
    }

    return ({
        'TRIMMER': [module.tool for module in tool_chain_set if isinstance(module, Trimmer)],
        'FASTQCER': [module.tool for module in tool_chain_set if isinstance(module, Fastqcer)],
        'ALIGNER': [module.tool for module in tool_chain_set if isinstance(module, Aligner)],
        'DEDUPER': [module.tool for module in tool_chain_set if isinstance(module, DeDuper)],
        'CALIBRATOR': [module.tool for module in tool_chain_set if isinstance(module, Calibrator)],
        'COUNTER': [module.tool for module in tool_chain_set if isinstance(module, Counter)],
        'BAMStatist': (BAMStatist().tool_list if any(isinstance(module, BAMStatist)
                                                     for module in tool_chain_set) else [])
    })


def dev_fq2bedgraph_tool_ls(df: pd.DataFrame) -> dict[str, list[str] | list[LiteralString]]:
    tool_chain_set: set[PipelineModule] = {
        tool
        for _, row in df.iterrows()  # type: ignore
        for tool in Fq2BedGraphRun(row=row, root_dir=Path('result')).ModuleChain
    }

    return ({
        'TRIMMER': [module.tool for module in tool_chain_set if isinstance(module, Trimmer)],
        'FASTQCER': [module.tool for module in tool_chain_set if isinstance(module, Fastqcer)],
        'ALIGNER': [module.tool for module in tool_chain_set if isinstance(module, Aligner)],
        'DEDUPER': [module.tool for module in tool_chain_set if isinstance(module, DeDuper)],
        'CALIBRATOR': [module.tool for module in tool_chain_set if isinstance(module, Calibrator)],
        'COUNTER': [module.tool for module in tool_chain_set if isinstance(module, Counter)],
        'BAMSTATIST': (BAMStatist().tool_list if any(isinstance(module, BAMStatist)
                                                     for module in tool_chain_set) else [])
    })


# TODO: 添加bedgraph2dm相关解析，输出编码后的文件和文件夹名
class BedGraph2DMRun:
    DM_OUTPUT_MAP: dict[str, list[str]] = {
        'methylsig': ['bedgraph.zst']
    }

    def __init__(self, row: pd.Series, root_dir: Path):  # type: ignore
        for col in ('FASTQCER1', 'FASTQCER2', 'STATS1', 'STATS2'):
            row[col] = None

        self.A: Fq2BedGraphRun
        self.B: Fq2BedGraphRun

        (self.A,  # type: ignore
         self.B) = [Fq2BedGraphRun(pd.Series({col.replace(suffix, ''): row[col]  # type: ignore
                                              for col in row.columns if col.endswith(suffix)}),  # type: ignore
                                   root_dir=Path('./result'))
                    for suffix in ('1', '2')]
        self.DMer: str = str(row['DMER'])  # type: ignore
        self.Files: list[Path] = self._generate_output(root_dir=root_dir)  # type: ignore

    def _generate_output(self, root_dir: Path) -> list[Path]:
        # ! the work to generate input files for DM analysis is done in rule/{DMer}/prepare.smk
        return [root_dir / self.DMer / f'{self.A.BaseName}.{self.B.BaseName}' / f'{self.A.Code}.{self.B.Code}.{ext}'
                for ext in self.DM_OUTPUT_MAP[self.DMer]]


def read_bedgraph2dm_sheet(csv_path: Path,
                           pipeline_coding: dict = PIPELINE_CODING) -> pd.DataFrame:  # type: ignore
    query_string = ('FNAME1.notna() and TRIMMER1.notna() and '
                    'FNAME2.notna() and TRIMMER2.notna() and '
                    '(ALIGNER1.isna() | ALIGNER1.isin(@pipeline_coding["ALIGNER"].keys())) and '
                    '(ALIGNER2.isna() | ALIGNER2.isin(@pipeline_coding["ALIGNER"].keys())) and '
                    '(DEDUPER1.isna() | DEDUPER1.isin(@pipeline_coding["DEDUPER"].keys())) and '
                    '(DEDUPER2.isna() | DEDUPER2.isin(@pipeline_coding["DEDUPER"].keys())) and '
                    '(CALIBRATOR1.isna() | CALIBRATOR1.isin(@pipeline_coding["CALIBRATOR"].keys())) and '
                    '(CALIBRATOR2.isna() | CALIBRATOR2.isin(@pipeline_coding["CALIBRATOR"].keys())) and '
                    '(COUNTER1.isna() | COUNTER1.isin(@pipeline_coding["COUNTER"].keys())) and '
                    '(COUNTER2.isna() | COUNTER2.isin(@pipeline_coding["COUNTER"].keys())) and '
                    'DMER.isin(@pipeline_coding["DMER"].keys())')

    b2d_df: pd.DataFrame
    b2d_df = (pd.read_csv(csv_path, comment='#',  # type: ignore
                          dtype={'FNAME1': 'str', 'TRIMMER1': 'str', 'ALIGNER1': 'str',
                                 'DEDUPER1': 'str', 'CALIBRATOR1': 'str', 'COUNTER1': 'str',
                                 'FNAME2': 'str', 'TRIMMER2': 'str', 'ALIGNER2': 'str',
                                 'DEDUPER2': 'str', 'CALIBRATOR2': 'str', 'COUNTER2': 'str',
                                 'DMER': 'str'})
              .dropna(subset=['FNAME1', 'FNAME2', 'DMER'])
              .query(query_string, engine='python'))

    if b2d_df.empty:
        raise ValueError('No valid rows found in the sheet.')

    return b2d_df.replace(np.nan, None)  # type: ignore


def bedgraph2dm_file_ls(csv_path_str: str) -> list[str]:
    df = read_bedgraph2dm_sheet(csv_path=Path(csv_path_str))
    return ([
        file.relative_to(Path()).as_posix()
        for _, row in df.iterrows()  # type: ignore
        for file in BedGraph2DMRun(row=row, root_dir=Path.cwd() / 'result').Files
    ])


def bedgraph2dm_tool_ls(csv_path_str: str) -> set[str]:
    df = read_bedgraph2dm_sheet(csv_path=Path(csv_path_str))
    return ({
        BedGraph2DMRun(row=row, root_dir=Path.cwd() / 'result').DMer for _, row in df.iterrows()  # type: ignore
    })

# TODO: 添加新smk文件，解码文件夹和文件名，并运行差异甲基化分析


def main():
    arg_parser: ArgumentParser = ArgumentParser()
    arg_parser.add_argument('-i', '--input', type=str, required=True,
                            help='csv file to investigate')
    args: Namespace = arg_parser.parse_args()

    f2b_df: pd.DataFrame = read_fq2bedgraph_sheet(csv_path=Path(args.input))

    print(f2b_df)

    print(dev_fq2bedgraph_tool_ls(df=f2b_df))

    print(dev_fq2bedgraph_file_ls(df=f2b_df))


if __name__ == '__main__':
    main()
