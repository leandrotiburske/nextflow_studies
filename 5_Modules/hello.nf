#!/usr/bin/env nextflow

include { SPLITLETTERS ; CONVERTTOUPPER   } from './modules.nf'

params.greeting = 'Hello world!' 
greeting_ch = Channel.of(params.greeting) 

workflow  {
    greeting_ch = Channel.of(params.greeting)
    SPLITLETTERS(greeting_ch)
    CONVERTTOUPPER(SPLITLETTERS.out.flatten())
    CONVERTTOUPPER.out.view()
}